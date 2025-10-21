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
import yaml
from pathlib import Path
from datetime import datetime

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
    
    return flatten_dict(config)

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
    # Record start time for runtime tracking
    start_time = time.time()
    
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
  
  # Complete pipeline with guide selection
  python run_CRISPR_target_automation.py targets.txt --scan-pam --qc --select-guides
  
  # Note: Upstream/downstream distances and genome assembly are controlled by config.yaml
  
        """
    )
    
    # Required argument
    parser.add_argument("targets", help="Input file with coordinates")
    
    # UCSC sequence fetching options (controlled by config.yaml for reproducibility)
    parser.add_argument("--retries", type=int, default=int(CONFIG.get("UCSC_RETRIES", "3")), 
                       help="Retry attempts for UCSC requests (default: 3)")
    
    # PAM scanning and QC options
    parser.add_argument("--scan-pam", action="store_true", 
                       help="Scan sequences for SpCas9 PAM sites (NGG)")
    parser.add_argument("--qc", action="store_true", 
                       help="Apply quality control filters to PAM sites")
    # QC parameters are controlled by policy.yaml only for reproducibility
    
    # Output options
    # Output file names are now configured in config.yaml (no CLI overrides for manifest integrity)
    
    # Control options
    parser.add_argument("--select-guides", action="store_true", 
                       help="Select top guides based on policy parameters (default: disabled)")
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
    # Validate required config keys
    required_config_keys = ['UCSC_UPSTREAM_DISTANCE', 'UCSC_DOWNSTREAM_DISTANCE', 'UCSC_GENOME_ASSEMBLY']
    required_output_keys = ['OUTPUTS_UPSTREAM_SEQUENCES', 'OUTPUTS_DOWNSTREAM_SEQUENCES', 'OUTPUTS_CRISPR_CANDIDATES', 
                           'OUTPUTS_CRISPR_CANDIDATES_QC', 'OUTPUTS_TOP_GUIDES']
    missing_keys = [key for key in required_config_keys + required_output_keys if key not in CONFIG]
    if missing_keys:
        print(f"❌ Error: Missing required configuration keys in config.yaml:")
        for key in missing_keys:
            print(f"   - {key}")
        print(f"\n💡 Please add these keys to your config.yaml file.")
        sys.exit(1)
    
    print(f"Upstream distance: {CONFIG.get('UCSC_UPSTREAM_DISTANCE')} bp (from config.yaml)")
    print(f"Downstream distance: {CONFIG.get('UCSC_DOWNSTREAM_DISTANCE')} bp (from config.yaml)")
    print(f"Genome: {CONFIG.get('UCSC_GENOME_ASSEMBLY')} (from config.yaml)")
    if args.scan_pam:
        print("🔍 PAM scanning: ENABLED (SpCas9 NGG)")
    if args.qc:
        # Validate required QC policy keys
        required_qc_keys = ['POLICY_QUALITY_CONTROL_GC_MIN', 'POLICY_QUALITY_CONTROL_GC_MAX', 
                           'POLICY_QUALITY_CONTROL_MAX_POLY_T', 'POLICY_QUALITY_CONTROL_MAX_HOMOPOLYMER',
                           'POLICY_QUALITY_CONTROL_RESTRICTION_SITES']
        missing_qc_keys = [key for key in required_qc_keys if key not in CONFIG]
        if missing_qc_keys:
            print(f"❌ Error: Missing required QC policy keys in policy.yaml:")
            for key in missing_qc_keys:
                print(f"   - {key}")
            print(f"\n💡 Please add these keys to your policy.yaml file.")
            sys.exit(1)
        
        gc_min = CONFIG.get("POLICY_QUALITY_CONTROL_GC_MIN")
        gc_max = CONFIG.get("POLICY_QUALITY_CONTROL_GC_MAX")
        print(f"🔬 Quality control: ENABLED (GC {gc_min:.2f}-{gc_max:.2f} from policy.yaml)")
    print(f"📝 Log file: {log_file}")
    print()
    
    # Log pipeline configuration
    logger.info(f"Pipeline configuration:")
    logger.info(f"  Input file: {args.targets}")
    logger.info(f"  Upstream distance: {CONFIG.get('UCSC_UPSTREAM_DISTANCE')} bp")
    logger.info(f"  Downstream distance: {CONFIG.get('UCSC_DOWNSTREAM_DISTANCE')} bp")
    logger.info(f"  Genome: {CONFIG.get('UCSC_GENOME_ASSEMBLY')}")
    logger.info(f"  PAM scanning: {'ENABLED' if args.scan_pam else 'DISABLED'}")
    logger.info(f"  Quality control: {'ENABLED' if args.qc else 'DISABLED'}")
    if args.qc:
        gc_min = CONFIG.get("POLICY_QUALITY_CONTROL_GC_MIN")
        gc_max = CONFIG.get("POLICY_QUALITY_CONTROL_GC_MAX")
        max_poly_t = CONFIG.get("POLICY_QUALITY_CONTROL_MAX_POLY_T")
        max_homopolymer = CONFIG.get("POLICY_QUALITY_CONTROL_MAX_HOMOPOLYMER")
        logger.info(f"  QC parameters: GC {gc_min:.2f}-{gc_max:.2f}, max poly-T {max_poly_t}, max homopolymer {max_homopolymer}")
    
    # Step 1: Fetch sequences from UCSC with optional PAM scanning and QC
    ucsc_cmd = [
        "python", "utils/get_ucsc_sequences.py", args.targets,
        "--retries", str(args.retries)
    ]
    
    # Add PAM scanning and QC options if requested
    if args.scan_pam:
        ucsc_cmd.extend(["--scan-pam"])
    if args.qc:
        ucsc_cmd.extend(["--qc"])
    
    if not run_command(ucsc_cmd, "Fetching sequences from UCSC", logger):
        logger.error("UCSC sequence fetching failed. Pipeline terminated.")
        print("❌ UCSC sequence fetching failed. Exiting.")
        sys.exit(1)
    
    # Determine what files to analyze with IDT
    idt_files = []
    
    # Check for PAM candidates first (if PAM scanning was enabled)
    if args.scan_pam:
        candidates_file = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES")
        if Path(candidates_file).exists():
            logger.info(f"Found PAM candidates file: {candidates_file}")
            print(f"🔍 Found PAM candidates: {candidates_file}")
            idt_files.append(candidates_file)
        else:
            print(f"❌ PAM candidates file not found: {candidates_file}")
            sys.exit(1)
    # Otherwise check for sequence files
    else:
        upstream_file = CONFIG.get("OUTPUTS_UPSTREAM_SEQUENCES")
        downstream_file = CONFIG.get("OUTPUTS_DOWNSTREAM_SEQUENCES")
        
        upstream_exists = Path(upstream_file).exists()
        downstream_exists = Path(downstream_file).exists()
        
        logger.info(f"Checking sequence files: upstream={upstream_exists}, downstream={downstream_exists}")
        
        if upstream_exists:
            logger.info(f"Found upstream sequences: {upstream_file}")
            idt_files.append(upstream_file)
        if downstream_exists:
            logger.info(f"Found downstream sequences: {downstream_file}")
            idt_files.append(downstream_file)
    
    if not idt_files:
        error_msg = "No files found for IDT analysis"
        logger.error(error_msg)
        if args.scan_pam:
            candidates_file = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES")
            logger.error(f"Expected PAM candidates file: {candidates_file}")
            print("❌ Error: No files found for IDT analysis")
            print(f"   Expected: {candidates_file}")
        else:
            upstream_file = CONFIG.get("OUTPUTS_UPSTREAM_SEQUENCES")
            downstream_file = CONFIG.get("OUTPUTS_DOWNSTREAM_SEQUENCES")
            logger.error(f"Expected sequence files: {upstream_file} or {downstream_file}")
            print("❌ Error: No files found for IDT analysis")
            print(f"   Expected: {upstream_file} or {downstream_file}")
        sys.exit(1)
    
    logger.info(f"Files selected for IDT analysis: {idt_files}")
    
    # Step 2: Analyze with IDT
    idt_cmd = ["python", "utils/idt_batch_crispr.py"] + idt_files
    
    if not run_command(idt_cmd, "Analyzing sequences with IDT", logger):
        logger.error("IDT analysis failed. Pipeline terminated.")
        print("❌ IDT analysis failed. Check your session cookie in config.yaml")
        print(f"📝 Check log file for detailed error information: {log_file}")
        sys.exit(1)
    
    # Step 3: Select top guides based on policy (optional)
    if args.select_guides:
        # Validate required guide selection policy keys
        required_guide_keys = ['POLICY_GUIDE_SELECTION_MIN_ON_TARGET_SCORE', 'POLICY_GUIDE_SELECTION_MIN_OFF_TARGET_SCORE',
                              'POLICY_GUIDE_SELECTION_NUM_GUIDES_PER_GENE', 'POLICY_GUIDE_SELECTION_ACCEPTED_PAMS']
        missing_guide_keys = [key for key in required_guide_keys if key not in CONFIG]
        if missing_guide_keys:
            print(f"❌ Error: Missing required guide selection policy keys in policy.yaml:")
            for key in missing_guide_keys:
                print(f"   - {key}")
            print(f"\n💡 Please add these keys to your policy.yaml file.")
            sys.exit(1)
        
        idt_results = [f"{Path(f).stem}_idt.csv" for f in idt_files if Path(f"{Path(f).stem}_idt.csv").exists()]
        
        if idt_results:
            print(f"\n🎯 Selecting top guides from {len(idt_results)} IDT result files...")
            select_cmd = ["python", "utils/select_top_guides.py"] + idt_results
            
            if not run_command(select_cmd, "Selecting top guides", logger):
                logger.warning("Guide selection failed, but IDT results are still available")
                print("⚠️  Guide selection failed, but IDT results are still available")
            else:
                print("✅ Top guides selected successfully!")
        else:
            logger.warning("No IDT results found for guide selection")
            print("⚠️  No IDT results found for guide selection")
    else:
        logger.info("Guide selection skipped (use --select-guides to enable)")
        print("ℹ️  Guide selection skipped (use --select-guides to enable)")
    
    # Step 4: Summary
    print(f"\n{'='*60}")
    print("PIPELINE SUMMARY")
    print(f"{'='*60}")
    
    # Show PAM scanning results if enabled
    if args.scan_pam:
        candidates_file = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES")
        if Path(candidates_file).exists():
            with open(candidates_file, 'r') as f:
                pam_count = sum(1 for line in f if line.startswith('>'))
            print(f"🔍 PAM candidates: {candidates_file} ({pam_count} sites)")
        
        if args.qc:
            qc_file = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES_QC")
            if Path(qc_file).exists():
                with open(qc_file, 'r') as f:
                    qc_count = sum(1 for line in f) - 1  # Subtract header
                print(f"🔬 QC results: {qc_file} ({qc_count} candidates)")
    
    # Show IDT analysis results
    for idt_file in idt_files:
        idt_results = f"{Path(idt_file).stem}_idt.csv"
        if Path(idt_results).exists():
            with open(idt_results, 'r') as f:
                idt_count = sum(1 for line in f) - 1  # Subtract header
            print(f"✅ IDT results: {idt_results} ({idt_count} sequences)")
        else:
            print(f"⚠️  Input file: {idt_file} (no IDT results)")
    
    # Show top guides selection results
    top_guides_file = CONFIG.get("OUTPUTS_TOP_GUIDES")
    if Path(top_guides_file).exists():
        with open(top_guides_file, 'r') as f:
            top_count = sum(1 for line in f) - 1  # Subtract header
        print(f"🏆 Top guides: {top_guides_file} ({top_count} selected)")
    else:
        if not args.select_guides:
            print(f"ℹ️  Top guides: Not generated (use --select-guides to enable)")
        else:
            # Check if IDT results exist
            idt_results_exist = any(Path(f"{Path(f).stem}_idt.csv").exists() for f in idt_files)
            if not idt_results_exist:
                print(f"ℹ️  Top guides: Not generated (no IDT results available)")
            else:
                print(f"ℹ️  Top guides: Not generated (selection step failed)")
    
    # Cleanup intermediate files if requested
    if args.cleanup:
        print(f"\n🧹 Cleaning up intermediate files...")
        files_to_remove = [CONFIG.get("OUTPUTS_UPSTREAM_SEQUENCES"), CONFIG.get("OUTPUTS_DOWNSTREAM_SEQUENCES")]
        if args.scan_pam:
            files_to_remove.append(CONFIG.get("OUTPUTS_CRISPR_CANDIDATES"))
        for file in files_to_remove:
            if Path(file).exists():
                Path(file).unlink()
                print(f"   Removed: {file}")
    
    # Auto-generate manifest
    print(f"\n📋 Generating manifest for reproducibility...")
    try:
        from manifest import write_manifest
        import glob
        import pandas as pd
        
        # Auto-capture QC stats from specific QC files generated by this pipeline
        qc_files = []
        total_passed_qc = 0
        total_failed_qc = 0
        
        # Only check the QC file that was actually generated by this run
        if args.qc:
            qc_file = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES_QC")
            if Path(qc_file).exists():
                qc_files.append(qc_file)
        
        for qc_file in qc_files:
            try:
                df = pd.read_csv(qc_file)
                passed = len(df[df['qc_status'].str.startswith('Pass', na=False)])
                failed = len(df) - passed
                total_passed_qc += passed
                total_failed_qc += failed
            except Exception as e:
                print(f"   ⚠️  Could not read QC file {qc_file}: {e}")
        
        # Count PAM candidates if available
        pam_candidates_count = 0
        if args.scan_pam:
            candidates_file = CONFIG.get("OUTPUTS_CRISPR_CANDIDATES")
            if Path(candidates_file).exists():
                with open(candidates_file, 'r') as f:
                    pam_candidates_count = sum(1 for line in f if line.startswith('>'))
        
        # Count IDT results
        idt_results_count = 0
        for idt_file in idt_files:
            idt_results = f"{Path(idt_file).stem}_idt.csv"
            if Path(idt_results).exists():
                with open(idt_results, 'r') as f:
                    idt_results_count += sum(1 for line in f) - 1  # Subtract header
        
        # Calculate total runtime
        total_runtime_sec = round(time.time() - start_time, 2)
        
        # Get species mapping information for manifest
        from utils.idt_batch_crispr import species_map
        ucsc_assembly = CONFIG.get("UCSC_GENOME_ASSEMBLY")
        idt_species, idt_assembly = species_map.get(ucsc_assembly, ("unknown", "unknown"))
        
        # Generate manifest with all stats
        summary_stats = {
            "pipeline_type": "full_automation",
            "targets_processed": len(open(args.targets, 'r').readlines()),
            "pam_candidates_found": pam_candidates_count,
            "qc_files_found": len(qc_files),
            "total_passed_qc": total_passed_qc,
            "total_failed_qc": total_failed_qc,
            "qc_pass_rate": round(total_passed_qc/(total_passed_qc + total_failed_qc), 3) if (total_passed_qc + total_failed_qc) > 0 else 0,
            "idt_files_processed": len(idt_files),
            "idt_results_generated": idt_results_count,
            "genome_assembly": CONFIG.get("UCSC_GENOME_ASSEMBLY"),
            "idt_species": idt_species,
            "idt_genome_assembly": idt_assembly,
            "species_mapping": f"{ucsc_assembly} → {idt_species} ({idt_assembly})",
            "upstream_distance": CONFIG.get("UCSC_UPSTREAM_DISTANCE"),
            "downstream_distance": CONFIG.get("UCSC_DOWNSTREAM_DISTANCE"),
            "pam_scanning_enabled": args.scan_pam,
            "pam_pattern": CONFIG.get("PAM_SCANNING_PATTERN") if args.scan_pam else None,
            "qc_enabled": args.qc,
            "total_runtime_sec": total_runtime_sec
        }
        
        manifest_file = f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}_manifest.json"
        # Use absolute paths to ensure we find the config files
        config_path = Path(__file__).parent.parent / "config.yaml"
        policy_path = Path(__file__).parent.parent / "policy.yaml"
        
        # Import and call manifest generation
        from manifest import write_manifest
        write_manifest(config_path=str(config_path), policy_path=str(policy_path), 
                      summary_stats=summary_stats, output=manifest_file)
        
    except ImportError:
        print("   ℹ️  Install manifest.py for automatic run tracking")
    except Exception as e:
        print(f"   ⚠️  Could not generate manifest: {e}")
    
    logger.info("Pipeline completed successfully")
    print(f"\n🎉 Pipeline complete!")
    if args.scan_pam:
        print("💡 Tip: Check the QC results for quality-filtered candidates, then analyze the best ones with IDT.")
    print("💡 Tip: Open the *_idt.csv files in Excel and sort by 'on_plus_off' column to see best targets first.")
    print(f"📝 Detailed log saved to: {log_file}")

if __name__ == "__main__":
    main()