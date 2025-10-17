#!/usr/bin/env python3
"""
get_ucsc_sequences.py

Fetch upstream and downstream sequences for genomic coordinates
from the UCSC Genome Browser API. Optionally scan for SpCas9 PAM sites
and apply quality control filters.

Usage:
    python get_ucsc_sequences.py targets.txt
    python get_ucsc_sequences.py chr14:103928378-103928397:+ --scan-pam
    python get_ucsc_sequences.py chr17:7668402-7668421:+ --scan-pam --qc

Input file (targets.txt):
    Coordinate format only:
        chr14:103928378-103928397:+
        chr17:7668402-7668421:-
        chr1:1000000-1000020:+

Outputs:
    Upstream_sequences.txt
    Downstream_sequences.txt
    CRISPR_candidates.txt (if --scan-pam is used, FASTA format for IDT)
    CRISPR_candidates_qc.csv (if --scan-pam --qc is used, QC results)
"""

import argparse, requests, time, re, sys, yaml
from pathlib import Path
from pam_scanner import scan_spcas9_sites, write_crispr_candidates
from qc_ucsc_seq import basic_qc, qc_pam_sites

def load_config():
    """Load configuration from config.yaml and policy.yaml files."""
    config_file = Path(__file__).parent.parent / "config.yaml"
    if not config_file.exists():
        return {}
    
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

UCSC_BASE = "https://api.genome.ucsc.edu"


def fetch_sequence(chrom, start, end, strand="+", genome=None, max_retries=None):
    """Fetch DNA sequence for a region from UCSC with retry logic."""
    # Always use config values (parameters will be None from our pipeline)
    genome = CONFIG.get("UCSC_GENOME_ASSEMBLY")
    max_retries = int(CONFIG.get("UCSC_RETRIES", 3))
    
    url = f"{UCSC_BASE}/getData/sequence"
    params = {"genome": genome, "chrom": chrom, "start": start, "end": end}
    
    for attempt in range(max_retries):
        try:
            r = requests.get(url, params=params, timeout=30)
            r.raise_for_status()
            seq = r.json().get("dna", "")
            if strand == "-":
                comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
                seq = seq.translate(comp)[::-1]
            return seq.upper()
        except Exception as e:
            if attempt < max_retries - 1:
                print(f"⚠️  Attempt {attempt + 1} failed for {chrom}:{start}-{end}: {e}")
                print(f"    Retrying in {2 ** attempt} seconds...")
                time.sleep(2 ** attempt)  # Exponential backoff
            else:
                print(f"❌  Failed to fetch {chrom}:{start}-{end} after {max_retries} attempts: {e}")
                return ""

def parse_target(line):
    """Parse coordinate (chr:start-end[:strand]) or return None."""
    line = line.strip()
    m = re.match(r"^(chr[\w]+):(\d+)-(\d+)(?::([+-]))?$", line)
    if m:
        chrom, start, end, strand = m.groups()
        return chrom, int(start), int(end), strand or "+"
    return None


def get_flanking_sequences(chrom, start, end, up=None, down=None, strand="+", genome=None):
    """Fetch upstream and downstream sequences."""
    # Always use config values (parameters will be None from our pipeline)
    up = int(CONFIG.get("UCSC_UPSTREAM_DISTANCE"))
    down = int(CONFIG.get("UCSC_DOWNSTREAM_DISTANCE"))
    genome = CONFIG.get("UCSC_GENOME_ASSEMBLY")
    
    upstream_start = max(0, start - up)
    downstream_end = end + down
    up_seq = fetch_sequence(chrom, upstream_start, start, strand, genome)
    down_seq = fetch_sequence(chrom, end, downstream_end, strand, genome)
    return up_seq, down_seq

def write_pam_fasta(candidates, filename, qc_results=None):
    """Write PAM candidates to FASTA format for IDT analysis."""
    if not candidates:
        print(f"⚠️  No PAM candidates written to {filename}")
        return
    
    # If QC results provided, only include candidates that passed QC
    if qc_results:
        passed_candidates = []
        qc_dict = {r[1]: r[5] for r in qc_results}  # name -> qc_status
        for candidate in candidates:
            parent, name, spacer, pam, strand = candidate
            if qc_dict.get(name, "").startswith("Pass"):
                passed_candidates.append((name, spacer))
        candidates_to_write = passed_candidates
        print(f"🔬 QC filtering: {len(passed_candidates)}/{len(candidates)} candidates passed QC")
    else:
        # Write all candidates
        candidates_to_write = [(name, spacer) for parent, name, spacer, pam, strand in candidates]
    
    with open(filename, "w") as f:
        for name, spacer in candidates_to_write:
            f.write(f">{name}\n{spacer}\n")
    print(f"✅ Saved {filename} ({len(candidates_to_write)} PAM sequences)")

def main():
    ap = argparse.ArgumentParser(description="Fetch upstream/downstream sequences for genomic coordinates from UCSC.")
    ap.add_argument("input", nargs='?', help="Input file with coordinates, or a single coordinate.")
    # Upstream/downstream distances are controlled by config.yaml only for reproducibility
    # Genome assembly is controlled by config.yaml only for reproducibility
    # Output file names are now configured in config.yaml (no CLI overrides for manifest integrity)
    ap.add_argument("--retries", type=int, default=int(CONFIG.get("UCSC_RETRIES", "3")), 
                   help="Number of retry attempts for failed requests.")
    ap.add_argument("--scan-pam", action="store_true", 
                   help="Scan sequences for SpCas9 PAM sites and output CRISPR candidates.")
    
    # QC options (parameters controlled by policy.yaml for reproducibility)
    ap.add_argument("--qc", action="store_true",
                   help="Apply quality control filters to CRISPR candidates.")
    args = ap.parse_args()

    # Handle input - either file or single coordinate
    if not args.input:
        print("Usage: python get_ucsc_sequences.py <input_file> [options]")
        print("   or: python get_ucsc_sequences.py <coordinate> [options]")
        print("Examples:")
        print("  python get_ucsc_sequences.py targets.txt")
        print("  python get_ucsc_sequences.py chr17:7668402-7668421:+")
        print("  python get_ucsc_sequences.py chr14:103928378-103928397:-")
        print("  python get_ucsc_sequences.py targets.txt --scan-pam --qc")
        print("")
        print("Note: Upstream/downstream distances and genome assembly are controlled by config.yaml")
        return
    
    # Check if input is a file or a single coordinate
    if Path(args.input).exists():
        # It's a file - use config output names
        targets = Path(args.input).read_text().strip().splitlines()
        upstream_out = CONFIG.get("OUTPUTS_UPSTREAM_SEQUENCES")
        downstream_out = CONFIG.get("OUTPUTS_DOWNSTREAM_SEQUENCES")
    else:
        # It's a single coordinate - use config output names
        targets = [args.input]
        upstream_out = CONFIG.get("OUTPUTS_UPSTREAM_SEQUENCES")
        downstream_out = CONFIG.get("OUTPUTS_DOWNSTREAM_SEQUENCES")
    upstream_records, downstream_records = [], []
    crispr_candidates = []
    qc_candidates = []

    for i, line in enumerate(targets, 1):
        if not line.strip():
            continue
        
        # Skip comment lines (starting with #)
        if line.strip().startswith('#'):
            continue

        parsed = parse_target(line)
        if parsed:
            chrom, start, end, strand = parsed
            label = line
        else:
            print(f"❌ Invalid coordinate format: {line}")
            print("Expected format: chr:start-end:strand (e.g., chr17:7668402-7668421:+)")
            continue

        print(f"🔎 Fetching {label} ...")
        
        # Validate required config keys
        required_keys = ['UCSC_GENOME_ASSEMBLY', 'UCSC_UPSTREAM_DISTANCE', 'UCSC_DOWNSTREAM_DISTANCE']
        missing_keys = [key for key in required_keys if key not in CONFIG]
        if missing_keys:
            print(f"❌ Error: Missing required configuration keys in config.yaml:")
            for key in missing_keys:
                print(f"   - {key}")
            print(f"\n💡 Please add these keys to your config.yaml file.")
            sys.exit(1)
        
        genome = CONFIG.get("UCSC_GENOME_ASSEMBLY")
        up_distance = int(CONFIG.get("UCSC_UPSTREAM_DISTANCE"))
        down_distance = int(CONFIG.get("UCSC_DOWNSTREAM_DISTANCE"))
        up_seq, down_seq = get_flanking_sequences(chrom, start, end, up_distance, down_distance, strand, genome)
        
        if up_seq:
            upstream_records.append((f"{label}_upstream", up_seq))
            # Scan upstream sequence for PAM sites if requested
            if args.scan_pam:
                sites = scan_spcas9_sites(up_seq)
                for idx, (spacer, pam, site_strand, pos) in enumerate(sites, 1):
                    # Sanitize label for IDT API (remove colons and plus signs, keep coordinate hyphens)
                    safe_label = label.replace(":", "_").replace("+", "plus")
                    candidate = (
                        f"{safe_label}_upstream",
                        f"{safe_label}_upstream_g{idx}",
                        spacer,
                        pam,
                        site_strand
                    )
                    crispr_candidates.append(candidate)
                    
                    # Add to QC candidates list if QC is requested
                    if args.qc:
                        qc_candidates.append(candidate)
        
        if down_seq:
            downstream_records.append((f"{label}_downstream", down_seq))
            # Scan downstream sequence for PAM sites if requested
            if args.scan_pam:
                sites = scan_spcas9_sites(down_seq)
                for idx, (spacer, pam, site_strand, pos) in enumerate(sites, 1):
                    # Sanitize label for IDT API (remove colons and plus signs, keep coordinate hyphens)
                    safe_label = label.replace(":", "_").replace("+", "plus")
                    candidate = (
                        f"{safe_label}_downstream",
                        f"{safe_label}_downstream_g{idx}",
                        spacer,
                        pam,
                        site_strand
                    )
                    crispr_candidates.append(candidate)
                    
                    # Add to QC candidates list if QC is requested
                    if args.qc:
                        qc_candidates.append(candidate)
        
        time.sleep(1.0)  # UCSC API rate limit

    def write_fasta(records, filename):
        if not records:
            print(f"⚠️  No sequences written to {filename}")
            return
        with open(filename, "w") as f:
            for name, seq in records:
                f.write(f">{name}\n{seq}\n")
        print(f"✅ Saved {filename} ({len(records)} sequences)")



    write_fasta(upstream_records, upstream_out)
    write_fasta(downstream_records, downstream_out)
    
    if args.scan_pam:
        if args.qc:
            # Apply QC to all candidates using policy parameters
            qc_results = qc_pam_sites(qc_candidates)
            
            # Write QC results (CSV format)
            qc_output = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES_QC")
            with open(qc_output, "w") as f:
                f.write("parent,name,spacer,pam,strand,qc_status\n")
                for result in qc_results:
                    parent, name, spacer, pam, strand, qc_status = result
                    f.write(f"{parent},{name},{spacer},{pam},{strand},{qc_status}\n")
            
            passed = sum(1 for r in qc_results if r[5].startswith("Pass"))
            total = len(qc_results)
            qc_output = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES_QC")
            print(f"✅ Saved {qc_output} ({passed}/{total} candidates passed QC)")
            
            # Write PAM candidates in FASTA format (QC-filtered)
            candidates_out = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES")
            write_pam_fasta(crispr_candidates, candidates_out, qc_results)
        else:
            # Write PAM candidates in FASTA format (all candidates)
            candidates_out = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES")
            write_pam_fasta(crispr_candidates, candidates_out)
    
    # Encourage manifest creation
    print("\n📋 Run completed! Consider creating a manifest for reproducibility:")
    print(f"   python manifest.py --config config.yaml --policy policy.yaml --stats '{{\"targets_processed\": {len(targets)}, \"pam_sites_found\": {len(crispr_candidates)}, \"genome_assembly\": \"{CONFIG.get('UCSC_GENOME_ASSEMBLY')}\"}}'")

if __name__ == "__main__":
    main()