#!/usr/bin/env python3
"""
manifest.py

Manifest generation for CRISPR pipeline runs.
Provides traceability and reproducibility by recording run metadata,
configuration hashes, and summary statistics.
"""

import json
import yaml
import hashlib
import datetime
import platform
import getpass
import socket
from pathlib import Path


def file_hash(path):
    """Generate MD5 hash of a file (first 7 characters for brevity)."""
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()[:7]


def write_manifest(config_path="../config.yaml", policy_path="../policy.yaml", summary_stats=None, output="manifest.json"):
    """
    Write a manifest file with run metadata and configuration hashes.
    
    Args:
        config_path: Path to config.yaml file (default: "config.yaml")
        policy_path: Path to policy.yaml file (default: "policy.yaml")
        summary_stats: Dict with summary statistics (optional)
        output: Output manifest file path (default: "manifest.json")
    """
    
    # Convert to Path objects for robust path resolution
    config_file = Path(config_path)
    policy_file = Path(policy_path)
    
    with open(config_file) as f: 
        cfg = yaml.safe_load(f)
    with open(policy_file) as f: 
        policy = yaml.safe_load(f)

    manifest = {
        "run_id": datetime.datetime.now(datetime.timezone.utc).isoformat().replace('+00:00', 'Z'),
        "guideforge_version": "1.1.0",
        "backend": "IDT",
        "genome_assembly": cfg.get("ucsc", {}).get("genome_assembly"),
        "config_hash": file_hash(config_file),
        "policy_hash": file_hash(policy_file),
        "python_version": platform.python_version(),
        "user": getpass.getuser(),
        "hostname": socket.gethostname()
    }
    
    # Add summary stats if provided
    if summary_stats:
        manifest.update(summary_stats)

    with open(output, "w") as out:
        json.dump(manifest, out, indent=2)
    
    print(f"📋 Manifest written to: {output}")
    print(f"   Run ID: {manifest['run_id']}")
    print(f"   Config hash: {manifest['config_hash']}")
    print(f"   Policy hash: {manifest['policy_hash']}")
    if summary_stats:
        print(f"   Summary stats: {len(summary_stats)} metrics recorded")


def read_manifest(manifest_path="manifest.json"):
    """Read and return manifest data from a JSON file."""
    with open(manifest_path, "r") as f:
        return json.load(f)


def compare_manifests(manifest1_path, manifest2_path):
    """
    Compare two manifest files and report differences.
    Useful for tracking changes between pipeline runs.
    """
    manifest1 = read_manifest(manifest1_path)
    manifest2 = read_manifest(manifest2_path)
    
    print("🔍 Manifest Comparison:")
    print(f"Run 1: {manifest1['run_id']} (User: {manifest1['user']})")
    print(f"Run 2: {manifest2['run_id']} (User: {manifest2['user']})")
    print()
    
    # Compare key fields
    key_fields = ["config_hash", "policy_hash", "python_version", "backend", "genome_assembly"]
    
    for field in key_fields:
        val1 = manifest1.get(field, "N/A")
        val2 = manifest2.get(field, "N/A")
        status = "✅" if val1 == val2 else "❌"
        print(f"{status} {field}: {val1} vs {val2}")
    
    # Compare summary stats
    print("\n📊 Summary Statistics:")
    all_keys = set(manifest1.keys()) | set(manifest2.keys())
    summary_keys = [k for k in all_keys if k.startswith(('n_', 'total_', 'passed_', 'failed_', 'pam_', 'qc_', 'idt_', 'pipeline_', 'genome_', 'upstream_', 'downstream_'))]
    
    for key in sorted(summary_keys):
        val1 = manifest1.get(key, "N/A")
        val2 = manifest2.get(key, "N/A")
        if val1 != val2:
            print(f"📈 {key}: {val1} → {val2}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="CRISPR Pipeline Manifest Tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create a basic manifest
  python manifest.py --config config.yaml --policy policy.yaml --output manifest.json
  
  # Create manifest with summary stats
  python manifest.py --config config.yaml --policy policy.yaml --stats '{"n_guides": 1000, "passed_qc": 850}'
  
  # Create manifest with custom policy file
  python manifest.py --config config.yaml --policy policy_strict.yaml --output strict_manifest.json
  
  # Compare two manifests
  python manifest.py --compare manifest1.json manifest2.json
        """
    )
    
    parser.add_argument("--config", default="config.yaml", 
                       help="Config file path (default: config.yaml)")
    parser.add_argument("--policy", default="policy.yaml", 
                       help="Policy file path (default: policy.yaml)")
    parser.add_argument("--output", default="manifest.json", 
                       help="Output manifest file (default: manifest.json)")
    parser.add_argument("--stats", 
                       help="Summary statistics as JSON string (e.g. '{\"n_guides\": 1000}')")
    parser.add_argument("--compare", nargs=2, metavar=("MANIFEST1", "MANIFEST2"), 
                       help="Compare two manifest files")
    
    args = parser.parse_args()
    
    if args.compare:
        compare_manifests(args.compare[0], args.compare[1])
    else:
        # Parse summary stats if provided
        summary_stats = None
        if args.stats:
            try:
                summary_stats = json.loads(args.stats)
            except json.JSONDecodeError as e:
                print(f"❌ Error parsing stats JSON: {e}")
                exit(1)
        
        write_manifest(args.config, args.policy, summary_stats, args.output)
