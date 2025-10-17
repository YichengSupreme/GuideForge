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

import argparse, re, sys, yaml
from pathlib import Path

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

# -------------------------------
# 🧪 Individual QC Checks
# -------------------------------

def gc_content(seq: str) -> float:
    """Return GC content as a fraction (0–1)."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) if len(seq) > 0 else 0.0


def has_poly_t(seq: str, length: int = None) -> bool:
    """Check if sequence contains a run of 'T's that may terminate transcription."""
    # Always use config values (parameters will be None from our pipeline)
    length = int(CONFIG.get("POLICY_QUALITY_CONTROL_MAX_POLY_T"))
    return "T" * length in seq.upper()


def has_homopolymer(seq: str, max_len: int = None) -> bool:
    """Detect any homopolymer (AAAAA, CCCCC, etc.) longer than allowed."""
    # Always use config values (parameters will be None from our pipeline)
    max_len = int(CONFIG.get("POLICY_QUALITY_CONTROL_MAX_HOMOPOLYMER"))
    seq = seq.upper()
    return bool(re.search(r"(A{%d,}|T{%d,}|C{%d,}|G{%d,})" % (max_len, max_len, max_len, max_len), seq))


def has_restriction_site(seq: str) -> bool:
    """
    Detect restriction sites that may interfere with plasmid assembly.
    Sites are configurable via policy.yaml.
    """
    seq = seq.upper()
    
    # Get restriction sites from policy config (no defaults - must be defined in policy.yaml)
    if "POLICY_QUALITY_CONTROL_RESTRICTION_SITES" not in CONFIG:
        print(f"❌ Error: Missing required policy key 'POLICY_QUALITY_CONTROL_RESTRICTION_SITES' in policy.yaml")
        print(f"💡 Please add this key to your policy.yaml file.")
        sys.exit(1)
    
    restriction_sites = CONFIG.get("POLICY_QUALITY_CONTROL_RESTRICTION_SITES")
    
    return any(site in seq for site in restriction_sites)


def has_excluded_motifs(seq: str) -> bool:
    """
    Check if sequence contains excluded motifs (e.g., poly-T, poly-A, poly-G).
    Motifs are configurable via policy.yaml.
    """
    seq = seq.upper()
    
    # Get excluded motifs from policy config
    excluded_motifs = CONFIG.get("POLICY_FILTERS_EXCLUDE_MOTIFS", [])
    
    return any(motif in seq for motif in excluded_motifs)


def gc_within_range(seq: str, gc_min: float = None, gc_max: float = None) -> bool:
    """Check if GC content is within an acceptable range for stable binding."""
    # Use policy values if not provided
    if gc_min is None:
        gc_min = float(CONFIG.get("POLICY_QUALITY_CONTROL_GC_MIN"))
    if gc_max is None:
        gc_max = float(CONFIG.get("POLICY_QUALITY_CONTROL_GC_MAX"))
    
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
    
    # Use policy config values as defaults if not provided
    required_qc_keys = ['POLICY_QUALITY_CONTROL_GC_MIN', 'POLICY_QUALITY_CONTROL_GC_MAX', 
                       'POLICY_QUALITY_CONTROL_MAX_POLY_T', 'POLICY_QUALITY_CONTROL_MAX_HOMOPOLYMER']
    missing_keys = [key for key in required_qc_keys if key not in CONFIG]
    if missing_keys:
        print(f"❌ Error: Missing required QC policy keys in policy.yaml:")
        for key in missing_keys:
            print(f"   - {key}")
        print(f"💡 Please add these keys to your policy.yaml file.")
        sys.exit(1)
    
    gc_min = gc_min or float(CONFIG.get("POLICY_QUALITY_CONTROL_GC_MIN"))
    gc_max = gc_max or float(CONFIG.get("POLICY_QUALITY_CONTROL_GC_MAX"))
    max_poly_t = max_poly_t or int(CONFIG.get("POLICY_QUALITY_CONTROL_MAX_POLY_T"))
    max_homopolymer = max_homopolymer or int(CONFIG.get("POLICY_QUALITY_CONTROL_MAX_HOMOPOLYMER"))

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

    # Excluded motifs
    if has_excluded_motifs(seq):
        return False, "Excluded motif"

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
    # QC parameters are now controlled by policy.yaml only for reproducibility
    # Use --policy flag to specify different policy files if needed
    ap.add_argument("--filtered-only", action="store_true", help="Only output candidates that pass QC")
    
    args = ap.parse_args()
    
    if not Path(args.input).exists():
        print(f"❌ Input file not found: {args.input}")
        return
    
    # Get QC parameters from policy (no defaults - must be defined in policy.yaml)
    gc_min = float(CONFIG.get("POLICY_QUALITY_CONTROL_GC_MIN"))
    gc_max = float(CONFIG.get("POLICY_QUALITY_CONTROL_GC_MAX"))
    max_poly_t = int(CONFIG.get("POLICY_QUALITY_CONTROL_MAX_POLY_T"))
    max_homopolymer = int(CONFIG.get("POLICY_QUALITY_CONTROL_MAX_HOMOPOLYMER"))
    
    print(f"🔬 QC Parameters (from policy.yaml): GC {gc_min:.2f}-{gc_max:.2f}, max poly-T {max_poly_t}, max homopolymer {max_homopolymer}")
    
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
    
    # Apply QC using policy parameters
    qc_results = qc_pam_sites(
        candidates,
        gc_min=gc_min,
        gc_max=gc_max,
        max_poly_t=max_poly_t,
        max_homopolymer=max_homopolymer
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