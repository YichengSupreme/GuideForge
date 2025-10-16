#!/usr/bin/env python3
"""
qc_ucsc_seq.py

Simple quality control (QC) filters for CRISPR gRNA sequences.

Functions:
    basic_qc(seq): Check a single sequence
    qc_pam_sites(candidates): QC a list of PAM sites
    main(): Standalone CLI for QC'ing PAM sites from file

Usage:
    # As module
    from qc_ucsc_seq import basic_qc, qc_pam_sites
    is_valid, reason = basic_qc("ACACTCATTGCAGACTCAGG")
    
    # Standalone
    python qc_ucsc_seq.py CRISPR_candidates.txt --output qc_results.txt
"""

import argparse, re, sys
from pathlib import Path

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

# -------------------------------
# 🧪 Individual QC Checks
# -------------------------------

def gc_content(seq: str) -> float:
    """Return GC content as a fraction (0–1)."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) if len(seq) > 0 else 0.0


def has_poly_t(seq: str, length: int = 4) -> bool:
    """Check if sequence contains a run of 'T's that may terminate transcription."""
    return "T" * length in seq.upper()


def has_homopolymer(seq: str, max_len: int = 5) -> bool:
    """Detect any homopolymer (AAAAA, CCCCC, etc.) longer than allowed."""
    seq = seq.upper()
    return bool(re.search(r"(A{%d,}|T{%d,}|C{%d,}|G{%d,})" % (max_len, max_len, max_len, max_len), seq))


def has_restriction_site(seq: str) -> bool:
    """
    Detect restriction sites that may interfere with plasmid assembly.
    Sites are configurable via config.sh QC_RESTRICTION_SITES setting.
    """
    seq = seq.upper()
    
    # Get restriction sites from config, fallback to defaults
    restriction_sites_str = CONFIG.get("QC_RESTRICTION_SITES", "GAATTC,AAGCTT,GGATCC,GGTACC,GCGGCCGC")
    restriction_sites = [site.strip() for site in restriction_sites_str.split(",") if site.strip()]
    
    return any(site in seq for site in restriction_sites)


def gc_within_range(seq: str, gc_min: float = 0.35, gc_max: float = 0.80) -> bool:
    """Check if GC content is within an acceptable range for stable binding."""
    gc = gc_content(seq)
    return gc_min <= gc <= gc_max


# -------------------------------
# 🧩 Combined QC
# -------------------------------

def basic_qc(seq: str,
             gc_min: float = None,
             gc_max: float = None,
             max_poly_t: int = None,
             max_homopolymer: int = None) -> tuple[bool, str]:
    """
    Run a basic QC pipeline on a single gRNA sequence.

    Returns:
        (is_valid, reason)
        is_valid: True if sequence passes all filters
        reason:   "Pass" or description of the first failing check
    """
    seq = seq.upper()
    
    # Use config values as defaults if not provided
    gc_min = gc_min or float(CONFIG.get("QC_GC_MIN", "0.35"))
    gc_max = gc_max or float(CONFIG.get("QC_GC_MAX", "0.80"))
    max_poly_t = max_poly_t or int(CONFIG.get("QC_MAX_POLY_T", "4"))
    max_homopolymer = max_homopolymer or int(CONFIG.get("QC_MAX_HOMOPOLYMER", "5"))

    # GC content check
    gc = gc_content(seq)
    if gc < gc_min:
        return False, f"Low GC ({gc:.2f})"
    if gc > gc_max:
        return False, f"High GC ({gc:.2f})"

    # Transcription terminator
    if has_poly_t(seq, max_poly_t):
        return False, f"PolyT (>{max_poly_t})"

    # Homopolymers (general)
    if has_homopolymer(seq, max_homopolymer):
        return False, f"Homopolymer (>{max_homopolymer})"

    # Restriction enzyme sites
    if has_restriction_site(seq):
        return False, "Restriction site"

    return True, "Pass"


# -------------------------------
# 🧩 Main QC Function
# -------------------------------

def qc_pam_sites(candidates, gc_min=None, gc_max=None, max_poly_t=None, max_homopolymer=None):
    """Apply QC filters to a list of PAM sites.
    
    Args:
        candidates (list): List of tuples (parent, name, spacer, pam, strand)
        gc_min (float): Minimum GC content
        gc_max (float): Maximum GC content
        max_poly_t (int): Maximum consecutive T's
        max_homopolymer (int): Maximum homopolymer length
        
    Returns:
        list: List of tuples (parent, name, spacer, pam, strand, qc_status)
    """
    qc_results = []
    
    for candidate in candidates:
        parent, name, spacer, pam, strand = candidate
        is_valid, reason = basic_qc(
            spacer,
            gc_min=gc_min,
            gc_max=gc_max,
            max_poly_t=max_poly_t,
            max_homopolymer=max_homopolymer
        )
        
        qc_status = reason if is_valid else f"FAIL: {reason}"
        qc_results.append((parent, name, spacer, pam, strand, qc_status))
    
    return qc_results

def main():
    """Simple standalone QC for PAM sites."""
    ap = argparse.ArgumentParser(description="Apply QC filters to CRISPR PAM sites")
    ap.add_argument("input", help="Input CSV file with PAM sites")
    ap.add_argument("--output", default="qc_results.txt", help="Output file (default: qc_results.txt)")
    ap.add_argument("--gc-min", type=float, default=float(CONFIG.get("QC_GC_MIN", "0.35")), 
                   help=f"Minimum GC content (default: {CONFIG.get('QC_GC_MIN', '0.35')})")
    ap.add_argument("--gc-max", type=float, default=float(CONFIG.get("QC_GC_MAX", "0.80")), 
                   help=f"Maximum GC content (default: {CONFIG.get('QC_GC_MAX', '0.80')})")
    ap.add_argument("--max-poly-t", type=int, default=int(CONFIG.get("QC_MAX_POLY_T", "4")), 
                   help=f"Maximum consecutive T's (default: {CONFIG.get('QC_MAX_POLY_T', '4')})")
    ap.add_argument("--max-homopolymer", type=int, default=int(CONFIG.get("QC_MAX_HOMOPOLYMER", "5")), 
                   help=f"Maximum homopolymer length (default: {CONFIG.get('QC_MAX_HOMOPOLYMER', '5')})")
    ap.add_argument("--filtered-only", action="store_true", help="Only output candidates that pass QC")
    
    args = ap.parse_args()
    
    if not Path(args.input).exists():
        print(f"❌ Input file not found: {args.input}")
        return
    
    print(f"🔬 QC Parameters: GC {args.gc_min:.2f}-{args.gc_max:.2f}")
    
    # Read PAM sites
    candidates = []
    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('parent'):  # Skip header
                parts = line.split(',')
                if len(parts) >= 5:
                    parent, name, spacer, pam, strand = parts[:5]
                    candidates.append((parent, name, spacer, pam, strand))
    
    print(f"📖 Found {len(candidates)} PAM sites")
    
    # Apply QC
    qc_results = qc_pam_sites(
        candidates,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        max_poly_t=args.max_poly_t,
        max_homopolymer=args.max_homopolymer
    )
    
    # Filter if requested
    if args.filtered_only:
        qc_results = [r for r in qc_results if r[5].startswith("Pass")]
        print(f"   Filtered to {len(qc_results)} candidates that passed QC")
    
    # Write results
    with open(args.output, "w") as f:
        f.write("parent,name,spacer,pam,strand,qc_status\n")
        for result in qc_results:
            parent, name, spacer, pam, strand, qc_status = result
            f.write(f"{parent},{name},{spacer},{pam},{strand},{qc_status}\n")
    
    passed = sum(1 for r in qc_results if r[5].startswith("Pass"))
    print(f"✅ Saved {args.output} ({passed}/{len(qc_results)} candidates passed QC)")


# -------------------------------
# 🧪 Example (for testing)
# -------------------------------

if __name__ == "__main__":
    main()