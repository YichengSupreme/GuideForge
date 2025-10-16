#!/usr/bin/env python3
"""
CRISPR Target Automation Pipeline

A complete pipeline that:
1. Fetches sequences from UCSC Genome Browser
2. Scans for SpCas9 PAM sites (NGG)
3. Applies quality control filters
4. Analyzes sequences with IDT CRISPR tools
5. Generates ranked results for target selection

Usage:
    python run_CRISPR_target_automation.py targets.txt [options]

Examples:
    python run_CRISPR_target_automation.py targets.txt
    python run_CRISPR_target_automation.py targets.txt --up 200 --down 200 --qc
    python run_CRISPR_target_automation.py targets.txt --genome hg19 --scan-pam --qc
"""

import argparse
import subprocess
import sys
import time
import logging
from pathlib import Path
from datetime import datetime

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

def setup_logging(log_file=None):
    """Set up logging to both file and console."""
    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = f"crispr_pipeline_{timestamp}.log"
    
    # Create logger
    logger = logging.getLogger('crispr_pipeline')
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

def run_command(cmd, description, logger):
    """Run a command and handle errors."""
    logger.info(f"Starting step: {description}")
    logger.debug(f"Command: {' '.join(cmd)}")
    
    print(f"\n{'='*60}")
    print(f"STEP: {description}")
    print(f"{'='*60}")
    print(f"Running: {' '.join(cmd)}")
    print()
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Log successful execution
        logger.info(f"Step '{description}' completed successfully")
        logger.debug(f"Command output:\n{result.stdout}")
        
        print(result.stdout)
        if result.stderr:
            logger.warning(f"Command stderr: {result.stderr}")
            print("Warnings/Info:", result.stderr)
        
        return True
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Step '{description}' failed with return code {e.returncode}"
        logger.error(error_msg)
        logger.error(f"Command: {' '.join(cmd)}")
        logger.error(f"Return code: {e.returncode}")
        logger.error(f"Error output: {e.stderr}")
        logger.error(f"Standard output: {e.stdout}")
        
        print(f"❌ Error in {description}:")
        print(f"Command: {' '.join(cmd)}")
        print(f"Return code: {e.returncode}")
        print(f"Error output: {e.stderr}")
        print(f"Standard output: {e.stdout}")
        
        return False
        
    except Exception as e:
        error_msg = f"Unexpected error in step '{description}': {str(e)}"
        logger.error(error_msg, exc_info=True)
        print(f"❌ Unexpected error in {description}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Complete CRISPR target automation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (sequences only)
  python run_CRISPR_target_automation.py targets.txt
  
  # With PAM scanning
  python run_CRISPR_target_automation.py targets.txt --scan-pam
  
  # Full pipeline with PAM scanning and QC
  python run_CRISPR_target_automation.py targets.txt --scan-pam --qc
  
  # Custom parameters
  python run_CRISPR_target_automation.py targets.txt --up 200 --down 200 --scan-pam --qc --gc-min 0.4
  
  # Different genome
  python run_CRISPR_target_automation.py targets.txt --genome hg19 --scan-pam --qc
  
        """
    )
    
    # Required argument
    parser.add_argument("targets", help="Input file with coordinates")
    
    # UCSC sequence fetching options
    parser.add_argument("--up", type=int, default=int(CONFIG.get("UPSTREAM_DISTANCE", "100")), 
                       help="Upstream distance in bp (default: 100)")
    parser.add_argument("--down", type=int, default=int(CONFIG.get("DOWNSTREAM_DISTANCE", "100")), 
                       help="Downstream distance in bp (default: 100)")
    parser.add_argument("--genome", default=CONFIG.get("GENOME_ASSEMBLY", "hg38"), 
                       help="Genome assembly (default: hg38)")
    parser.add_argument("--retries", type=int, default=int(CONFIG.get("UCSC_RETRIES", "3")), 
                       help="Retry attempts for UCSC requests (default: 3)")
    
    # PAM scanning and QC options
    parser.add_argument("--scan-pam", action="store_true", 
                       help="Scan sequences for SpCas9 PAM sites (NGG)")
    parser.add_argument("--qc", action="store_true", 
                       help="Apply quality control filters to PAM sites")
    parser.add_argument("--gc-min", type=float, default=0.35, 
                       help="Minimum GC content for QC (default: 0.35)")
    parser.add_argument("--gc-max", type=float, default=0.80, 
                       help="Maximum GC content for QC (default: 0.80)")
    parser.add_argument("--max-poly-t", type=int, default=4, 
                       help="Maximum consecutive T's for QC (default: 4)")
    parser.add_argument("--max-homopolymer", type=int, default=5, 
                       help="Maximum homopolymer length for QC (default: 5)")
    
    # Output options
    parser.add_argument("--upstream-out", default=CONFIG.get("UPSTREAM_OUTPUT", "Upstream_sequences.txt"), 
                       help="Upstream output filename (default: Upstream_sequences.txt)")
    parser.add_argument("--downstream-out", default=CONFIG.get("DOWNSTREAM_OUTPUT", "Downstream_sequences.txt"), 
                       help="Downstream output filename (default: Downstream_sequences.txt)")
    parser.add_argument("--candidates-out", default="CRISPR_candidates.txt",
                       help="PAM candidates output filename (default: CRISPR_candidates.txt)")
    parser.add_argument("--qc-output", default="CRISPR_candidates_qc.csv",
                       help="QC results output filename (default: CRISPR_candidates_qc.csv)")
    
    # Control options
    parser.add_argument("--cleanup", action="store_true", 
                       help="Remove intermediate files after completion")
    parser.add_argument("--log-file", 
                       help="Specify log file name (default: auto-generated with timestamp)")
    
    args = parser.parse_args()
    
    # Set up logging
    logger, log_file = setup_logging(args.log_file)
    logger.info("Starting CRISPR Target Automation Pipeline")
    logger.info(f"Command line arguments: {sys.argv}")
    logger.info(f"Configuration loaded: {len(CONFIG)} settings")
    
    # Validate input file
    if not Path(args.targets).exists():
        error_msg = f"Input file '{args.targets}' not found"
        logger.error(error_msg)
        print(f"❌ Error: {error_msg}")
        sys.exit(1)
    
    logger.info(f"Input file validated: {args.targets}")
    
    print("🧬 CRISPR Target Automation Pipeline")
    print("=" * 50)
    print(f"Input file: {args.targets}")
    print(f"Upstream distance: {args.up} bp")
    print(f"Downstream distance: {args.down} bp")
    print(f"Genome: {args.genome}")
    if args.scan_pam:
        print("🔍 PAM scanning: ENABLED (SpCas9 NGG)")
    if args.qc:
        print(f"🔬 Quality control: ENABLED (GC {args.gc_min:.2f}-{args.gc_max:.2f})")
    print(f"📝 Log file: {log_file}")
    print()
    
    # Log pipeline configuration
    logger.info(f"Pipeline configuration:")
    logger.info(f"  Input file: {args.targets}")
    logger.info(f"  Upstream distance: {args.up} bp")
    logger.info(f"  Downstream distance: {args.down} bp")
    logger.info(f"  Genome: {args.genome}")
    logger.info(f"  PAM scanning: {'ENABLED' if args.scan_pam else 'DISABLED'}")
    logger.info(f"  Quality control: {'ENABLED' if args.qc else 'DISABLED'}")
    if args.qc:
        logger.info(f"  QC parameters: GC {args.gc_min:.2f}-{args.gc_max:.2f}, max poly-T {args.max_poly_t}, max homopolymer {args.max_homopolymer}")
    
    # Step 1: Fetch sequences from UCSC with optional PAM scanning and QC
    ucsc_cmd = [
        "python", "get_ucsc_sequences.py", args.targets,
        "--up", str(args.up),
        "--down", str(args.down),
        "--genome", args.genome,
        "--retries", str(args.retries),
        "--upstream-out", args.upstream_out,
        "--downstream-out", args.downstream_out
    ]
    
    # Add PAM scanning and QC options if requested
    if args.scan_pam:
        ucsc_cmd.extend(["--scan-pam", "--candidates-out", args.candidates_out])
    if args.qc:
        ucsc_cmd.extend([
            "--qc",
            "--gc-min", str(args.gc_min),
            "--gc-max", str(args.gc_max),
            "--max-poly-t", str(args.max_poly_t),
            "--max-homopolymer", str(args.max_homopolymer),
            "--qc-output", args.qc_output
        ])
    
    if not run_command(ucsc_cmd, "Fetching sequences from UCSC", logger):
        logger.error("UCSC sequence fetching failed. Pipeline terminated.")
        print("❌ UCSC sequence fetching failed. Exiting.")
        sys.exit(1)
    
    # Determine what files to analyze with IDT
    idt_files = []
    
    # Check for PAM candidates first (if PAM scanning was enabled)
    if args.scan_pam and Path(args.candidates_out).exists():
        logger.info(f"Found PAM candidates file: {args.candidates_out}")
        print(f"🔍 Found PAM candidates: {args.candidates_out}")
        idt_files.append(args.candidates_out)
    # Otherwise check for sequence files
    else:
        upstream_exists = Path(args.upstream_out).exists()
        downstream_exists = Path(args.downstream_out).exists()
        
        logger.info(f"Checking sequence files: upstream={upstream_exists}, downstream={downstream_exists}")
        
        if upstream_exists:
            logger.info(f"Found upstream sequences: {args.upstream_out}")
            idt_files.append(args.upstream_out)
        if downstream_exists:
            logger.info(f"Found downstream sequences: {args.downstream_out}")
            idt_files.append(args.downstream_out)
    
    if not idt_files:
        error_msg = "No files found for IDT analysis"
        logger.error(error_msg)
        if args.scan_pam:
            logger.error(f"Expected PAM candidates file: {args.candidates_out}")
            print("❌ Error: No files found for IDT analysis")
            print(f"   Expected: {args.candidates_out}")
        else:
            logger.error(f"Expected sequence files: {args.upstream_out} or {args.downstream_out}")
            print("❌ Error: No files found for IDT analysis")
            print(f"   Expected: {args.upstream_out} or {args.downstream_out}")
        sys.exit(1)
    
    logger.info(f"Files selected for IDT analysis: {idt_files}")
    
    # Step 2: Analyze with IDT
    idt_cmd = ["python", "idt_batch_crispr.py"] + idt_files
    
    if not run_command(idt_cmd, "Analyzing sequences with IDT", logger):
        logger.error("IDT analysis failed. Pipeline terminated.")
        print("❌ IDT analysis failed. Check your session cookie in config.sh")
        print(f"📝 Check log file for detailed error information: {log_file}")
        sys.exit(1)
    
    # Step 3: Summary
    print(f"\n{'='*60}")
    print("PIPELINE SUMMARY")
    print(f"{'='*60}")
    
    # Show PAM scanning results if enabled
    if args.scan_pam:
        if Path(args.candidates_out).exists():
            with open(args.candidates_out, 'r') as f:
                pam_count = sum(1 for line in f if line.startswith('>'))
            print(f"🔍 PAM candidates: {args.candidates_out} ({pam_count} sites)")
        
        if args.qc and Path(args.qc_output).exists():
            with open(args.qc_output, 'r') as f:
                qc_count = sum(1 for line in f) - 1  # Subtract header
            print(f"🔬 QC results: {args.qc_output} ({qc_count} candidates)")
    
    # Show IDT analysis results
    for idt_file in idt_files:
        idt_results = f"{Path(idt_file).stem}_idt.csv"
        if Path(idt_results).exists():
            with open(idt_results, 'r') as f:
                idt_count = sum(1 for line in f) - 1  # Subtract header
            print(f"✅ IDT results: {idt_results} ({idt_count} sequences)")
        else:
            print(f"⚠️  Input file: {idt_file} (no IDT results)")
    
    # Cleanup intermediate files if requested
    if args.cleanup:
        print(f"\n🧹 Cleaning up intermediate files...")
        files_to_remove = [args.upstream_out, args.downstream_out]
        if args.scan_pam:
            files_to_remove.append(args.candidates_out)
        for file in files_to_remove:
            if Path(file).exists():
                Path(file).unlink()
                print(f"   Removed: {file}")
    
    logger.info("Pipeline completed successfully")
    print(f"\n🎉 Pipeline complete!")
    if args.scan_pam:
        print("💡 Tip: Check the QC results for quality-filtered candidates, then analyze the best ones with IDT.")
    print("💡 Tip: Open the *_idt.csv files in Excel and sort by 'on_plus_off' column to see best targets first.")
    print(f"📝 Detailed log saved to: {log_file}")

if __name__ == "__main__":
    main()
