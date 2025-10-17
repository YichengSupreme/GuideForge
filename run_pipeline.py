#!/usr/bin/env python3
"""
CRISPR Target Automation Pipeline Launcher

This is a simple launcher script that runs the main pipeline from the utils directory.
All scripts have been moved to the utils/ folder for better organization.

Usage:
    python run_pipeline.py targets.txt [options]
    
Examples:
    python run_pipeline.py targets.txt
    python run_pipeline.py targets.txt --scan-pam --qc
    python run_pipeline.py targets.txt --scan-pam --qc --select-guides
"""

import sys
import subprocess
from pathlib import Path

def main():
    """Launch the main pipeline script from the utils directory."""
    # Get the path to the utils directory
    utils_dir = Path(__file__).parent / "utils"
    pipeline_script = utils_dir / "run_CRISPR_target_automation.py"
    
    # Check if the pipeline script exists
    if not pipeline_script.exists():
        print("❌ Error: Pipeline script not found in utils/ directory")
        print(f"   Expected: {pipeline_script}")
        sys.exit(1)
    
    # Run the pipeline script with all arguments passed through
    cmd = ["python", str(pipeline_script)] + sys.argv[1:]
    
    try:
        result = subprocess.run(cmd, check=True)
        sys.exit(result.returncode)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)
    except KeyboardInterrupt:
        print("\n⚠️  Pipeline interrupted by user")
        sys.exit(1)

if __name__ == "__main__":
    main()
