#!/bin/bash
# CRISPR Pipeline Configuration
# Edit these values as needed for your lab setup

# UCSC Sequence Fetching Settings (get_ucsc_sequences.py)
UPSTREAM_DISTANCE=100
DOWNSTREAM_DISTANCE=100
GENOME_ASSEMBLY="hg38"
UCSC_RETRIES=3
UPSTREAM_OUTPUT="Upstream_sequences.txt"
DOWNSTREAM_OUTPUT="Downstream_sequences.txt"

# IDT Analysis Settings (idt_batch_crispr.py)
IDT_BATCH_SIZE=10
IDT_TIMEOUT=60
IDT_SESSION_COOKIE=""
# Output Settings
UPSTREAM_RESULTS="Upstream_sequences_idt.csv"
DOWNSTREAM_RESULTS="Downstream_sequences_idt.csv"

# Quality Control Settings (qc_ucsc_seq.py)
QC_GC_MIN=0.35
QC_GC_MAX=0.80
QC_MAX_POLY_T=4
QC_MAX_HOMOPOLYMER=5
# Restriction sites to avoid (comma-separated): EcoRI, HindIII, BamHI, KpnI, NotI
QC_RESTRICTION_SITES="GAATTC,AAGCTT,GGATCC,GGTACC,GCGGCCGC"

# System Settings 
PYTHON_CMD="python3"
VERBOSE=1
