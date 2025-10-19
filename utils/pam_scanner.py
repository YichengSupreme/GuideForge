#!/usr/bin/env python3
"""
pam_scanner.py

PAM site scanning functions for CRISPR target identification.
Scans DNA sequences for SpCas9 PAM sites (NGG) and returns spacer sequences.

Functions:
    reverse_complement(seq): Generate reverse complement of DNA sequence
    scan_spcas9_sites(seq): Scan sequence for SpCas9 PAM sites
    scan_sequences(sequences, names): Scan multiple sequences with names
    write_crispr_candidates(candidates, filename): Write candidates to CSV

Usage:
    from pam_scanner import scan_spcas9_sites, scan_sequences
    sites = scan_spcas9_sites("ATCGATCGATCGATCGATCGNGG")
"""

import re
from pathlib import Path


def reverse_complement(seq):
    """Generate reverse complement of DNA sequence.
    
    Args:
        seq (str): DNA sequence
        
    Returns:
        str: Reverse complement sequence
    """
    tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tbl)[::-1]


def scan_spcas9_sites(seq, pam_pattern):
    """Scan sequence for PAM sites with configurable PAM pattern.
    
    Args:
        seq (str): DNA sequence to scan
        pam_pattern (str): PAM pattern (e.g., "NGG", "RGG", "VGG")
        
    Returns:
        list: List of tuples (spacer, pam, strand, position)
    """
    # Define nucleotide ambiguity codes
    ambiguity_codes = {
        'N': '[ACGT]',      # Any nucleotide
        'R': '[AG]',        # Purine (A or G)
        'Y': '[CT]',        # Pyrimidine (C or T)
        'S': '[GC]',        # Strong (G or C)
        'W': '[AT]',        # Weak (A or T)
        'K': '[GT]',        # Keto (G or T)
        'M': '[AC]',        # Amino (A or C)
        'B': '[CGT]',       # Not A
        'D': '[AGT]',       # Not C
        'H': '[ACT]',       # Not G
        'V': '[ACG]',       # Not T
    }
    
    # Convert PAM pattern to regex
    pam_regex = ""
    for char in pam_pattern:
        if char in ambiguity_codes:
            pam_regex += ambiguity_codes[char]
        else:
            pam_regex += char
    
    sites = []
    
    # + strand: 20nt spacer followed by PAM
    for m in re.finditer(f'(?=([ACGT]{{20}})({pam_regex}))', seq):
        spacer, pam = m.group(1), m.group(2)
        sites.append((spacer, pam, '+', m.start()))
    
    # - strand: scan reverse complement, then report spacer in genomic (+) orientation
    rc_seq = reverse_complement(seq)
    seq_len = len(seq)
    for m in re.finditer(f'(?=([ACGT]{{20}})({pam_regex}))', rc_seq):
        spacer_rc, pam_rc = m.group(1), m.group(2)
        spacer = reverse_complement(spacer_rc)  # report 5'->3' on genomic +
        pam = reverse_complement(pam_rc)
        # Convert position back to genomic coordinates
        pos = seq_len - (m.start() + 20 + len(pam_pattern))
        sites.append((spacer, pam, '-', pos))
    
    return sites


def scan_sequences(sequences, names, pam_pattern):
    """Scan multiple sequences for PAM sites.
    
    Args:
        sequences (list): List of DNA sequences
        names (list): List of sequence names
        pam_pattern (str): PAM pattern (e.g., "NGG", "RGG", "VGG")
        
    Returns:
        list: List of tuples (parent, name, spacer, pam, strand)
    """
    candidates = []
    
    for seq, name in zip(sequences, names):
        sites = scan_spcas9_sites(seq, pam_pattern)
        for idx, (spacer, pam, strand, pos) in enumerate(sites, 1):
            candidates.append((
                name,
                f"{name}_g{idx}",
                spacer,
                pam,
                strand
            ))
    
    return candidates


def write_crispr_candidates(candidates, filename):
    """Write CRISPR candidates to CSV file.
    
    Args:
        candidates (list): List of candidate tuples
        filename (str): Output filename
    """
    if not candidates:
        print(f"⚠️  No CRISPR candidates found")
        return
    
    with open(filename, "w") as f:
        f.write("parent,name,spacer,pam,strand\n")
        for parent, name, spacer, pam, strand in candidates:
            f.write(f"{parent},{name},{spacer},{pam},{strand}\n")
    
    print(f"✅ Saved {filename} ({len(candidates)} candidates)")


def main():
    """Command-line interface for PAM scanning."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Scan DNA sequences for SpCas9 PAM sites")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("-o", "--output", default="CRISPR_candidates.txt",
                       help="Output CSV file (default: CRISPR_candidates.txt)")
    parser.add_argument("-p", "--pam-pattern", required=True,
                       help="PAM pattern (required). Supports ambiguity codes: N=any, R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC, B=CGT, D=AGT, H=ACT, V=ACG")
    args = parser.parse_args()
    
    # Read FASTA file
    sequences, names = [], []
    with open(args.input) as f:
        current_name, current_seq = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    sequences.append("".join(current_seq).upper())
                    names.append(current_name)
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            sequences.append("".join(current_seq).upper())
            names.append(current_name)
    
    # Scan for PAM sites
    candidates = scan_sequences(sequences, names, args.pam_pattern)
    
    # Write results
    write_crispr_candidates(candidates, args.output)


if __name__ == "__main__":
    main()
