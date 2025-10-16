#!/usr/bin/env python3
"""
get_ucsc_sequences.py

Fetch upstream and downstream sequences for genomic coordinates
from the UCSC Genome Browser API. Optionally scan for SpCas9 PAM sites
and apply quality control filters.

Usage:
    python get_ucsc_sequences.py targets.txt --up 100 --down 100 --genome hg38
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

import argparse, requests, time, re, sys
from pathlib import Path
from pam_scanner import scan_spcas9_sites, write_crispr_candidates
from qc_ucsc_seq import basic_qc, qc_pam_sites

def load_config():
    """Load configuration from config.sh file."""
    config_file = Path(__file__).parent / "config.sh"
    if not config_file.exists():
        return {}
    
    config = {}
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#') and '=' in line:
                key, value = line.split('=', 1)
                value = value.strip('"\'')
                config[key.strip()] = value
    return config

CONFIG = load_config()

UCSC_BASE = "https://api.genome.ucsc.edu"


def fetch_sequence(chrom, start, end, strand="+", genome="hg38", max_retries=3):
    """Fetch DNA sequence for a region from UCSC with retry logic."""
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


def get_flanking_sequences(chrom, start, end, up=100, down=100, strand="+", genome="hg38"):
    """Fetch upstream and downstream sequences."""
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
    ap.add_argument("--up", type=int, default=int(CONFIG.get("UPSTREAM_DISTANCE", "100")), 
                   help="Number of bp upstream to fetch.")
    ap.add_argument("--down", type=int, default=int(CONFIG.get("DOWNSTREAM_DISTANCE", "100")), 
                   help="Number of bp downstream to fetch.")
    ap.add_argument("--genome", default=CONFIG.get("GENOME_ASSEMBLY", "hg38"), 
                   help="Genome assembly (default: hg38).")
    ap.add_argument("--upstream-out", default=CONFIG.get("UPSTREAM_OUTPUT", "Upstream_sequences.txt"), 
                   help="Output filename for upstream sequences.")
    ap.add_argument("--downstream-out", default=CONFIG.get("DOWNSTREAM_OUTPUT", "Downstream_sequences.txt"), 
                   help="Output filename for downstream sequences.")
    ap.add_argument("--retries", type=int, default=int(CONFIG.get("UCSC_RETRIES", "3")), 
                   help="Number of retry attempts for failed requests.")
    ap.add_argument("--scan-pam", action="store_true", 
                   help="Scan sequences for SpCas9 PAM sites and output CRISPR candidates.")
    ap.add_argument("--candidates-out", default="CRISPR_candidates.txt",
                   help="Output filename for CRISPR candidates (default: CRISPR_candidates.txt)")
    
    # QC options
    ap.add_argument("--qc", action="store_true",
                   help="Apply quality control filters to CRISPR candidates.")
    ap.add_argument("--gc-min", type=float, default=0.35, 
                   help="Minimum GC content for QC (default: 0.35)")
    ap.add_argument("--gc-max", type=float, default=0.80, 
                   help="Maximum GC content for QC (default: 0.80)")
    ap.add_argument("--max-poly-t", type=int, default=4, 
                   help="Maximum consecutive T's for QC (default: 4)")
    ap.add_argument("--max-homopolymer", type=int, default=5, 
                   help="Maximum homopolymer length for QC (default: 5)")
    ap.add_argument("--qc-output", default="CRISPR_candidates_qc.csv",
                   help="Output filename for QC results (default: CRISPR_candidates_qc.csv)")
    args = ap.parse_args()

    # Handle input - either file or single coordinate
    if not args.input:
        print("Usage: python get_ucsc_sequences.py <input_file> [options]")
        print("   or: python get_ucsc_sequences.py <coordinate> [options]")
        print("Examples:")
        print("  python get_ucsc_sequences.py targets.txt")
        print("  python get_ucsc_sequences.py chr17:7668402-7668421:+")
        print("  python get_ucsc_sequences.py chr14:103928378-103928397:-")
        return
    
    # Check if input is a file or a single coordinate
    if Path(args.input).exists():
        # It's a file
        targets = Path(args.input).read_text().strip().splitlines()
        # Use default output names for file input
        upstream_out = args.upstream_out
        downstream_out = args.downstream_out
    else:
        # It's a single coordinate
        targets = [args.input]
        # Create descriptive filenames for single input
        safe_name = args.input.replace(":", "_").replace("-", "_").replace("+", "plus").replace("-", "minus")
        upstream_out = f"{safe_name}_upstream.txt"
        downstream_out = f"{safe_name}_downstream.txt"
    upstream_records, downstream_records = [], []
    crispr_candidates = []
    qc_candidates = []

    for i, line in enumerate(targets, 1):
        if not line.strip():
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
        up_seq, down_seq = get_flanking_sequences(chrom, start, end, args.up, args.down, strand, args.genome)
        
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
            # Apply QC to all candidates
            qc_results = qc_pam_sites(
                qc_candidates,
                gc_min=args.gc_min,
                gc_max=args.gc_max,
                max_poly_t=args.max_poly_t,
                max_homopolymer=args.max_homopolymer
            )
            
            # Write QC results (CSV format)
            with open(args.qc_output, "w") as f:
                f.write("parent,name,spacer,pam,strand,qc_status\n")
                for result in qc_results:
                    parent, name, spacer, pam, strand, qc_status = result
                    f.write(f"{parent},{name},{spacer},{pam},{strand},{qc_status}\n")
            
            passed = sum(1 for r in qc_results if r[5].startswith("Pass"))
            total = len(qc_results)
            print(f"✅ Saved {args.qc_output} ({passed}/{total} candidates passed QC)")
            
            # Write PAM candidates in FASTA format (QC-filtered)
            write_pam_fasta(crispr_candidates, args.candidates_out, qc_results)
        else:
            # Write PAM candidates in FASTA format (all candidates)
            write_pam_fasta(crispr_candidates, args.candidates_out)

if __name__ == "__main__":
    main()