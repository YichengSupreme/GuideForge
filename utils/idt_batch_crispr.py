#!/usr/bin/env python3
"""
IDT CRISPR batch analyzer (v2)

- Reads all .txt files in current directory (FASTA-style)
- Submits sequences to IDT's CRISPR /search/ API in batches
- Polls the /getresult/ endpoint until full on/off scores are available
- Saves results to <filename>_idt.csv

Usage:
    python idt_batch_crispr.py                    # Run connectivity test only
    python idt_batch_crispr.py file1.txt [file2.txt ...]  # Analyze sequences
"""

import os, time, random, requests, pandas as pd, sys, re, logging, yaml
from pathlib import Path
from datetime import datetime

# === 1. Load configuration ===
def load_config():
    config_file = Path(__file__).parent.parent / "config.yaml"
    if not config_file.exists():
        print("Error: config.yaml not found.")
        sys.exit(1)
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Load policy file if specified
    policy_file = config.get('policy_file', './policy.yaml')
    policy_path = Path(__file__).parent.parent / policy_file
    
    if policy_path.exists():
        with open(policy_path, 'r') as f:
            policy = yaml.safe_load(f)
        
        # Merge policy into config
        config['policy'] = policy
    
    # Flatten the nested structure for backward compatibility
    def flatten_dict(d, prefix=''):
        flattened = {}
        for key, value in d.items():
            new_key = f'{prefix}_{key.upper()}' if prefix else key.upper()
            if isinstance(value, dict):
                flattened.update(flatten_dict(value, new_key))
            else:
                flattened[new_key] = value
        return flattened
    
    flattened = flatten_dict(config)
    
    return flattened

CONFIG = load_config()
COOKIE_STRING = CONFIG.get("IDT_SESSION_COOKIE", "").strip()
ENDPOINT_SEARCH = "https://eu.idtdna.com/sciservices/sherlock/crispr/search/"
ENDPOINT_RESULT = "https://eu.idtdna.com/sciservices/sherlock/crispr/getresult/"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json, text/javascript, */*; q=0.01",
    "User-Agent": "Mozilla/5.0 (compatible; IDT-CRISPR-BatchBot/2.0)",
    "Cookie": COOKIE_STRING
}
# Validate required IDT config keys
required_idt_keys = ['IDT_BATCH_SIZE', 'IDT_TIMEOUT', 'IDT_DELAY_MIN', 'IDT_DELAY_MAX', 'IDT_POLL_INTERVAL', 'IDT_POLL_MAX', 'IDT_RETRY_ATTEMPTS']
missing_keys = [key for key in required_idt_keys if key not in CONFIG]
if missing_keys:
    print(f"❌ Error: Missing required IDT configuration keys in config.yaml:")
    for key in missing_keys:
        print(f"   - {key}")
    print(f"\n💡 Please add these keys to your config.yaml file.")
    sys.exit(1)

BATCH_SIZE = int(CONFIG.get("IDT_BATCH_SIZE"))
DELAY_MIN = int(CONFIG.get("IDT_DELAY_MIN"))
DELAY_MAX = int(CONFIG.get("IDT_DELAY_MAX"))
TIMEOUT = int(CONFIG.get("IDT_TIMEOUT"))
POLL_INTERVAL = int(CONFIG.get("IDT_POLL_INTERVAL"))
POLL_MAX = int(CONFIG.get("IDT_POLL_MAX"))
RETRY_ATTEMPTS = int(CONFIG.get("IDT_RETRY_ATTEMPTS"))

# === 1.5. Setup logging ===
def setup_logging():
    """Set up logging for IDT batch processing."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"idt_batch_{timestamp}.log"
    
    # Create logger
    logger = logging.getLogger('idt_batch')
    logger.setLevel(logging.DEBUG)
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # Create formatters
    file_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_formatter = logging.Formatter(
        '%(levelname)s - %(message)s'
    )
    
    # File handler
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    return logger, log_file

# === 2. FASTA parser ===
def read_fasta_txt(path, logger):
    """Read FASTA file and return list of (name, sequence) tuples."""
    logger.debug(f"Reading FASTA file: {path}")
    names, seqs = [], []
    
    try:
        with open(path) as f:
            current_name, current_seq = None, []
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_name:
                        seq = "".join(current_seq).upper()
                        if len(seq) == 0:
                            logger.warning(f"Empty sequence for {current_name} at line {line_num}")
                        else:
                            seqs.append(seq)
                            names.append(current_name)
                    current_name = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_name:
                seq = "".join(current_seq).upper()
                if len(seq) == 0:
                    logger.warning(f"Empty sequence for {current_name} at end of file")
                else:
                    seqs.append(seq)
                    names.append(current_name)
        
        logger.info(f"Successfully read {len(names)} sequences from {path}")
        return list(zip(names, seqs))
        
    except FileNotFoundError:
        logger.error(f"File not found: {path}")
        raise
    except Exception as e:
        logger.error(f"Error reading file {path}: {str(e)}")
        raise

# === 3a. Species / genome mapping ===
species_map = {
    "hg38": ("human", "hg38"),
    "hg19": ("human", "hg19"),
    "mm10": ("mouse", "mm10"),
    "mm39": ("mouse", "mm39"),
    "rn7": ("rat", "rn7"),
    "danRer11": ("zebrafish", "danRer11"),
    "ce11": ("c_elegans", "ce11"),
}

# Check genome assembly support first
ucsc_assembly = CONFIG.get("UCSC_GENOME_ASSEMBLY")
if not ucsc_assembly:
    print("❌ Error: UCSC_GENOME_ASSEMBLY not found in config.yaml")
    print("Please add 'genome_assembly' to the 'ucsc' section of your config.yaml")
    sys.exit(1)

if ucsc_assembly not in species_map:
    print(f"❌ Error: Unsupported genome assembly '{ucsc_assembly}'")
    print(f"Supported assemblies: {', '.join(species_map.keys())}")
    print(f"Please update your config.yaml with a supported assembly or add '{ucsc_assembly}' to the species_map")
    sys.exit(1)

idt_species, idt_assembly = species_map[ucsc_assembly]
print(f"✅ Genome mapping: UCSC {ucsc_assembly} → IDT {idt_species} ({idt_assembly})")

# === 3. Payload builder ===
def build_payload(named_seq_list, logger):
    """Build API payload for IDT CRISPR analysis."""
    logger.debug(f"Building payload for {len(named_seq_list)} sequences")
    
    # Validate sequences
    valid_sequences = []
    for name, seq in named_seq_list:
        if not seq or len(seq) < 10:
            logger.warning(f"Skipping sequence {name}: too short ({len(seq)} bp)")
            continue
        if not all(c in 'ACGTN' for c in seq.upper()):
            logger.warning(f"Skipping sequence {name}: contains invalid characters")
            continue
        valid_sequences.append((name, seq))
    
    if not valid_sequences:
        logger.error("No valid sequences found in batch")
        return None
    
    logger.debug(f"Using {len(valid_sequences)} valid sequences out of {len(named_seq_list)}")
    logger.debug(f"Using IDT species: {idt_species}, genome: {idt_assembly}")
    
    return {
        "ToolName": "CRISPR_SEQUENCE",
        "Types": [12],
        "DesignsPerGene": "6",
        "From": 1,
        "To": 1,
        "Species": idt_species,
        "Genome": idt_assembly,  # some IDT endpoints also accept Genome
        "NamedSequences": [{"Name": n, "Sequence": s} for n, s in valid_sequences]
    }

# === 4. Score helpers ===
def extract_numeric_score(score_text):
    if not score_text:
        return None
    nums = re.findall(r'-?\d+\.?\d*', str(score_text))
    return float(nums[0]) if nums else None

def calculate_on_off_score(on_val, off_val):
    on = extract_numeric_score(on_val)
    off = extract_numeric_score(off_val)
    if on is None or off is None:
        return None
    return on + off

# === 5. Parse utility ===
def parse_props(props):
    """Extract name, on-target, off-target from a 'Props' list."""
    name, on_val, off_val = None, None, None
    for p in props:
        if p.get("FieldName") == "SearchField":
            name = p.get("FieldValue")
        elif p.get("FieldName") == "CardMessages":
            for msg in p.get("FieldValue", []):
                dev, desc = msg.get("deviation"), msg.get("description")
                if dev == "BAD_FOR_POTENCY":
                    on_val = desc
                elif dev == "BAD_FOR_OFF_TARGET":
                    off_val = desc
        elif p.get("FieldName") in ["OnTargetPotential", "OffTargetRiskSpecificity"]:
            # These appear only in getresult/ endpoint
            if p.get("FieldName") == "OnTargetPotential":
                on_val = p.get("FieldValue")
            else:
                off_val = p.get("FieldValue")
    return name, on_val, off_val

def parse_response(resp_json, sequence_map=None, logger=None):
    """Parse IDT API response and extract results."""
    results = []
    try:
        top = resp_json[0] if isinstance(resp_json, list) else resp_json
        details = top.get("Details", [])
        logger.debug(f"Parsing response with {len(details)} details")
        
        for i, d in enumerate(details):
            name, on_val, off_val = parse_props(d.get("Props", []))
            if name:
                on_off = calculate_on_off_score(on_val, off_val)
                # Get the actual sequence if sequence_map is provided
                sequence = sequence_map.get(name, "") if sequence_map else ""
                results.append((name, sequence, on_val, off_val, on_off))
                logger.debug(f"Parsed result {i+1}: {name} - on={on_val}, off={off_val}, combined={on_off}")
            else:
                logger.warning(f"Detail {i+1} has no name, skipping")
                
    except Exception as e:
        error_msg = f"Parse error: {str(e)}"
        if logger:
            logger.error(error_msg, exc_info=True)
        else:
            print(f"⚠️ {error_msg}")
    return results

# === 6. Poll for completed results ===
def poll_for_result(lookup_key, sequence_map=None, logger=None):
    """Poll IDT getresult/<key> until scores appear."""
    url = ENDPOINT_RESULT + lookup_key
    logger.info(f"Polling for results: {url}")
    
    for attempt in range(POLL_MAX):
        try:
            logger.debug(f"Poll attempt {attempt + 1}/{POLL_MAX}")
            r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
            
            if r.status_code != 200:
                logger.warning(f"HTTP {r.status_code} on poll attempt {attempt + 1}")
                time.sleep(POLL_INTERVAL)
                continue
                
            data = r.json()
            parsed = parse_response(data, sequence_map, logger)
            
            # If we found actual numeric scores, return them
            if any(p[2] or p[3] for p in parsed):  # Updated indices for new column order
                logger.info(f"Found complete results after {attempt + 1} attempts")
                return parsed
            else:
                logger.debug(f"No complete scores yet, waiting {POLL_INTERVAL}s...")
                
            time.sleep(POLL_INTERVAL)
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Request error on poll attempt {attempt + 1}: {str(e)}")
            time.sleep(POLL_INTERVAL)
        except Exception as e:
            logger.error(f"Unexpected error on poll attempt {attempt + 1}: {str(e)}")
            time.sleep(POLL_INTERVAL)
    
    logger.warning(f"Polling timed out after {POLL_MAX} attempts")
    return []

# === 7. Tiny test ===
def tiny_test(logger):
    """Test IDT API connectivity with a single sequence."""
    logger.info("Running connectivity test with one sequence...")
    print("🔎 Running connectivity test with one sequence...")
    
    # Check if session cookie is configured
    cookie = CONFIG.get("IDT_SESSION_COOKIE", "").strip()
    if not cookie or cookie == "YOUR_IDT_SESSION_COOKIE_HERE":
        logger.error("IDT session cookie not configured")
        print("❌ Error: IDT session cookie not configured")
        print("Please update the 'session_cookie' in the 'idt' section of config.yaml")
        print("Get your session cookie from IDT's website after logging in")
        return False
    
    test_seq = "AACGCGCCGCGCGCCCTTGT"
    payload = build_payload([("TinyTest", test_seq)], logger)
    
    if not payload:
        logger.error("Failed to build test payload")
        print("❌ Failed to build test payload")
        return False
    
    try:
        r = requests.post(ENDPOINT_SEARCH, headers=HEADERS, json=payload, timeout=TIMEOUT)
        logger.info(f"Test request HTTP status: {r.status_code}")
        print("HTTP status:", r.status_code)
        
        if r.status_code == 401:
            logger.error("Authentication failed - invalid session cookie")
            print("❌ Authentication failed - invalid session cookie")
            print("Please check your session cookie in config.yaml")
            return False
        elif r.status_code == 500:
            logger.error(f"IDT server error: {r.text[:200]}")
            print("❌ IDT server error - this might be a temporary issue")
            print("Try again later or check IDT's service status")
            return False
        elif r.status_code != 200:
            logger.error(f"Connection failed with status {r.status_code}: {r.text[:200]}")
            print(f"❌ Connection failed with status {r.status_code}")
            print("Check your internet connection and IDT service status")
            return False
            
        data = r.json()
        lookup = data[0].get("LookupKey") if isinstance(data, list) else None
        sequence_map = {"TinyTest": test_seq}
        parsed = parse_response(data, sequence_map, logger)
        
        if not parsed and lookup:
            logger.info("No immediate results, polling for final result...")
            print("⌛ Polling for final result...")
            parsed = poll_for_result(lookup, sequence_map, logger)
        
        if parsed:
            logger.info(f"Test successful: {parsed}")
            print("✅ Test result:", parsed)
            return True
        else:
            logger.warning("Test completed but no results returned")
            print("⚠️ Test completed but no results returned")
            return False
            
    except Exception as e:
        logger.error(f"Test failed with error: {str(e)}", exc_info=True)
        print(f"❌ Test failed: {e}")
        return False

# === 8. Process file ===
def process_file(txt_path, logger):
    """Process a single FASTA file through IDT analysis."""
    name = Path(txt_path).stem
    logger.info(f"Processing file: {txt_path}")
    print(f"\n📄 Processing {name}...")
    
    try:
        named_seq = read_fasta_txt(txt_path, logger)
    except Exception as e:
        logger.error(f"Failed to read file {txt_path}: {str(e)}")
        print(f"❌ Failed to read file: {e}")
        return None
    
    # Create sequence map for lookup
    sequence_map = {name: seq for name, seq in named_seq}
    
    # Calculate total batches for progress tracking
    total_batches = (len(named_seq) + BATCH_SIZE - 1) // BATCH_SIZE
    total_sequences = len(named_seq)
    
    logger.info(f"Processing {total_sequences} sequences in {total_batches} batches")
    print(f"  Total sequences: {total_sequences} | Batches: {total_batches}")
    
    results = []
    for i in range(0, len(named_seq), BATCH_SIZE):
        batch_num = i // BATCH_SIZE + 1
        batch = named_seq[i:i+BATCH_SIZE]
        
        logger.info(f"Processing batch {batch_num}/{total_batches} with {len(batch)} sequences")
        print(f"  Batch {batch_num}/{total_batches}: Sending {len(batch)} sequences...")
        
        payload = build_payload(batch, logger)
        if not payload:
            logger.error(f"Failed to build payload for batch {batch_num}")
            print(f"    ❌ Failed to build payload for batch {batch_num}")
            continue
        
        # Retry logic for timeouts
        r = None
        for attempt in range(RETRY_ATTEMPTS):
            try:
                logger.debug(f"Batch {batch_num} attempt {attempt + 1}/{RETRY_ATTEMPTS}")
                r = requests.post(ENDPOINT_SEARCH, headers=HEADERS, json=payload, timeout=TIMEOUT)
                break
            except requests.exceptions.ReadTimeout:
                if attempt < RETRY_ATTEMPTS - 1:
                    logger.warning(f"Timeout on batch {batch_num} attempt {attempt + 1}, retrying...")
                    print(f"    ⏰ Timeout on attempt {attempt + 1}, retrying...")
                    time.sleep(5)
                else:
                    logger.error(f"Batch {batch_num} failed after {RETRY_ATTEMPTS} timeout attempts")
                    print(f"    ❌ Batch {batch_num} failed after {RETRY_ATTEMPTS} attempts")
                    break
            except Exception as e:
                logger.error(f"Error on batch {batch_num} attempt {attempt + 1}: {str(e)}")
                print(f"    ❌ Error on batch {batch_num}: {e}")
                break
        
        if r is None:
            continue
            
        if r.status_code != 200:
            logger.error(f"HTTP {r.status_code} for batch {batch_num}: {r.text[:300]}")
            print(f"    ❌ HTTP {r.status_code}: {r.text[:300]}")
            continue
            
        try:
            data = r.json()
            lookup = data[0].get("LookupKey") if isinstance(data, list) else None
            parsed = parse_response(data, sequence_map, logger)
            
            # If empty, poll for final results
            if (not parsed or all(p[2] is None for p in parsed)) and lookup:
                logger.info(f"Polling for final results for batch {batch_num}")
                print(f"    ⌛ Polling getresult/{lookup} ...")
                parsed = poll_for_result(lookup, sequence_map, logger)
            
            results.extend(parsed)
            logger.info(f"Batch {batch_num} completed with {len(parsed)} results")
            print(f"    ✅ Got {len(parsed)} results")
            
        except Exception as e:
            logger.error(f"Error processing response for batch {batch_num}: {str(e)}")
            print(f"    ❌ Error processing response: {e}")
            continue
        
        # Save progress every 10 batches
        if batch_num % 10 == 0:
            try:
                temp_df = pd.DataFrame(results, columns=["sequence_name", "dna_sequence", "on_target_score", "off_target_score", "on_plus_off"])
                temp_df = temp_df.sort_values("on_plus_off", ascending=False)
                temp_csv = f"{name}_idt_temp.csv"
                temp_df.to_csv(temp_csv, index=False)
                logger.info(f"Saved progress checkpoint: {temp_csv} ({len(results)} sequences)")
                print(f"    💾 Saved progress: {temp_csv} ({len(results)} sequences so far)")
            except Exception as e:
                logger.warning(f"Failed to save progress checkpoint: {str(e)}")
        
        time.sleep(random.uniform(DELAY_MIN, DELAY_MAX))
    
    # Save final results with grouping by parent sequence
    try:
        df = pd.DataFrame(results, columns=["sequence_name", "dna_sequence", "on_target_score", "off_target_score", "on_plus_off"])
        
        # Extract parent sequence from sequence_name (everything before the last underscore)
        df['parent_sequence'] = df['sequence_name'].str.rsplit('_', n=1).str[0]
        
        # Group by parent sequence and rank within each group by on_plus_off (descending)
        df['rank_within_parent'] = df.groupby('parent_sequence')['on_plus_off'].rank(method='dense', ascending=False).astype(int)
        
        # Sort by parent sequence, then by rank within parent
        df = df.sort_values(['parent_sequence', 'rank_within_parent'])
        
        # Reorder columns for better readability
        df = df[['parent_sequence', 'rank_within_parent', 'sequence_name', 'dna_sequence', 'on_target_score', 'off_target_score', 'on_plus_off']]
        
        out_csv = f"{name}_idt.csv"
        df.to_csv(out_csv, index=False)
        logger.info(f"Saved final results: {out_csv} ({len(df)} sequences)")
        print(f"✅ Saved {out_csv} ({len(df)} sequences)")
        
        # Print summary by parent sequence
        parent_summary = df.groupby('parent_sequence').size().sort_index()
        print(f"Results by gene:")
        for parent, count in parent_summary.items():
            print(f"   {parent}: {count} guides")
        
        # Clean up temp files if run was successful
        temp_csv = f"{name}_idt_temp.csv"
        if Path(temp_csv).exists():
            Path(temp_csv).unlink()
            logger.debug(f"Cleaned up temp file: {temp_csv}")
            print(f"🧹 Cleaned up {temp_csv}")
        
        return df
        
    except Exception as e:
        logger.error(f"Failed to save final results: {str(e)}")
        print(f"❌ Failed to save final results: {e}")
        return None

# === 9. Main ===
def main():
    # Set up logging first
    logger, log_file = setup_logging()
    logger.info("Starting IDT batch CRISPR analysis")
    logger.info(f"Command line arguments: {sys.argv}")
    
    if len(sys.argv) < 2:
        logger.info("No input files provided, running connectivity test only...")
        print("🔎 No input files provided, running connectivity test only...")
        print("📝 Log file:", log_file)
        
        # Run connectivity test
        if tiny_test(logger):
            logger.info("Connectivity test passed successfully")
            print("✅ Connectivity test passed! Your IDT session is working.")
            print("💡 To analyze sequences, run: python idt_batch_crispr.py file1.txt [file2.txt]")
        else:
            logger.error("Connectivity test failed")
            print("❌ Connectivity test failed. Check your session cookie in config.yaml")
            print(f"📝 Check log file for details: {log_file}")
        return
    
    input_files = [f for f in sys.argv[1:] if Path(f).suffix == ".txt" and Path(f).exists()]
    if not input_files:
        error_msg = "No valid input files found."
        logger.error(error_msg)
        print(f"❌ {error_msg}")
        return
    
    logger.info(f"Processing {len(input_files)} file(s): {', '.join(input_files)}")
    print(f"📁 Processing {len(input_files)} file(s): {', '.join(input_files)}")
    print(f"📝 Log file: {log_file}")
    
    # Run connectivity test
    if not tiny_test(logger):
        logger.error("Connectivity test failed. Exiting.")
        print("❌ Connectivity test failed. Check your session cookie in config.yaml")
        print(f"📝 Check log file for details: {log_file}")
        return
    
    # Process each file
    successful_files = 0
    for f in input_files:
        try:
            result = process_file(f, logger)
            if result is not None:
                successful_files += 1
                logger.info(f"Successfully processed {f}")
            else:
                logger.error(f"Failed to process {f}")
        except Exception as e:
            logger.error(f"Unexpected error processing {f}: {str(e)}", exc_info=True)
            print(f"❌ Unexpected error processing {f}: {e}")
    
    logger.info(f"Processing complete: {successful_files}/{len(input_files)} files successful")
    print(f"\n🎉 Processing complete: {successful_files}/{len(input_files)} files successful")
    print(f"📝 Detailed log saved to: {log_file}")
    
    # Encourage manifest creation
    print("\n📋 Run completed! Consider creating a manifest for reproducibility:")
    print(f"   python manifest.py --config config.yaml --policy policy.yaml --stats '{{\"files_processed\": {len(input_files)}, \"successful_files\": {successful_files}, \"success_rate\": {successful_files/len(input_files):.3f}, \"idt_species\": \"{idt_species}\", \"idt_genome_assembly\": \"{idt_assembly}\", \"species_mapping\": \"{ucsc_assembly} → {idt_species} ({idt_assembly})\"}}'")

if __name__ == "__main__":
    main()