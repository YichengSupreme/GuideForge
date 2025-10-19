#!/usr/bin/env python3
"""
select_top_guides.py

Selects top CRISPR guides from IDT analysis results based on policy parameters.
Applies guide selection criteria from policy.yaml to rank and filter candidates.

Usage:
    python select_top_guides.py [input_files...]
    
Examples:
    python select_top_guides.py CRISPR_candidates_idt.csv
    python select_top_guides.py file1_idt.csv file2_idt.csv
    python select_top_guides.py *_idt.csv
"""

import argparse
import pandas as pd
import yaml
import sys
from pathlib import Path
from datetime import datetime

def load_config():
    """Load configuration from config.yaml and policy.yaml files."""
    config_file = Path(__file__).parent.parent / "config.yaml"
    if not config_file.exists():
        print("Error: config.yaml not found.")
        sys.exit(1)
    
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

def load_idt_results(file_path):
    """Load IDT results from CSV file."""
    try:
        df = pd.read_csv(file_path)
        print(f"📖 Loaded {len(df)} sequences from {file_path}")
        return df
    except Exception as e:
        print(f"❌ Error loading {file_path}: {e}")
        return None

def apply_guide_selection_filters(df, config):
    """Apply guide selection filters based on policy parameters."""
    print("🔍 Applying guide selection filters...")
    
    original_count = len(df)
    
    # Get policy parameters
    min_on_target = config.get("POLICY_GUIDE_SELECTION_MIN_ON_TARGET_SCORE")
    min_off_target = config.get("POLICY_GUIDE_SELECTION_MIN_OFF_TARGET_SCORE")
    accepted_pams = config.get("POLICY_GUIDE_SELECTION_ACCEPTED_PAMS", [])
    
    # Filter by on-target score
    if min_on_target is not None:
        df = df[df['on_target_score'] >= min_on_target]
        print(f"   On-target score ≥ {min_on_target}: {len(df)} sequences")
    
    # Filter by off-target score (higher is better for IDT off_target_score)
    if min_off_target is not None:
        df = df[df['off_target_score'] >= min_off_target]
        print(f"   Off-target score ≥ {min_off_target}: {len(df)} sequences")
    
    # Filter by PAM sites (if PAM information is available)
    if accepted_pams and 'pam' in df.columns:
        df = df[df['pam'].isin(accepted_pams)]
        print(f"   Accepted PAMs {accepted_pams}: {len(df)} sequences")
    
    filtered_count = len(df)
    print(f"   Filtered: {original_count} → {filtered_count} sequences ({filtered_count/original_count*100:.1f}% passed)")
    
    return df

def select_top_guides(df, config):
    """Select top guides based on policy parameters."""
    print("🏆 Selecting top guides...")
    
    # Get policy parameters
    num_guides_per_gene = config.get("POLICY_GUIDE_SELECTION_NUM_GUIDES_PER_GENE", 5)
    min_spacing = config.get("POLICY_GUIDE_SELECTION_MIN_SPACING_BP", 30)
    
    # Sort by combined score (on_plus_off) descending
    df = df.sort_values('on_plus_off', ascending=False)
    
    # If we have parent sequence information, select top guides per parent
    if 'parent_sequence' in df.columns:
        top_guides = []
        
        for parent in df['parent_sequence'].unique():
            parent_guides = df[df['parent_sequence'] == parent].copy()
            
            # Select top guides for this parent sequence
            selected = []
            for _, guide in parent_guides.iterrows():
                if len(selected) >= num_guides_per_gene:
                    break
                
                # Check spacing from previously selected guides
                if min_spacing > 0 and len(selected) > 0:
                    # This is a simplified spacing check - in practice you'd need genomic coordinates
                    # For now, we'll just take the top N guides per parent
                    pass
                
                selected.append(guide)
            
            top_guides.extend(selected)
            print(f"   {parent}: selected {len(selected)} guides")
        
        result_df = pd.DataFrame(top_guides)
    else:
        # No parent information, just take top N overall
        result_df = df.head(num_guides_per_gene)
        print(f"   Selected top {len(result_df)} guides overall")
    
    return result_df

def main():
    parser = argparse.ArgumentParser(
        description="Select top CRISPR guides from IDT analysis results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single file
  python select_top_guides.py CRISPR_candidates_idt.csv
  
  # Process multiple files
  python select_top_guides.py file1_idt.csv file2_idt.csv
  
  # Process all IDT results
  python select_top_guides.py *_idt.csv
        """
    )
    
    parser.add_argument("input_files", nargs="+", help="IDT results CSV files to process")

    # Guide selection parameters are controlled by policy.yaml only for reproducibility
    
    args = parser.parse_args()
    
    print("🧬 CRISPR Guide Selection")
    print("=" * 50)
    
    # Load and combine all IDT results
    all_results = []
    for file_path in args.input_files:
        if not Path(file_path).exists():
            print(f"⚠️  File not found: {file_path}")
            continue
        
        df = load_idt_results(file_path)
        if df is not None:
            all_results.append(df)
    
    if not all_results:
        print("❌ No valid input files found")
        return
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    print(f"📊 Combined {len(combined_df)} sequences from {len(all_results)} files")
    
    # All guide selection parameters come from policy.yaml for reproducibility
    
    # Validate required guide selection policy keys (accepted_pams is optional)
    required_guide_keys = ['POLICY_GUIDE_SELECTION_MIN_ON_TARGET_SCORE', 'POLICY_GUIDE_SELECTION_MIN_OFF_TARGET_SCORE',
                          'POLICY_GUIDE_SELECTION_NUM_GUIDES_PER_GENE']
    missing_keys = [key for key in required_guide_keys if key not in CONFIG]
    if missing_keys:
        print(f"❌ Error: Missing required guide selection policy keys in policy.yaml:")
        for key in missing_keys:
            print(f"   - {key}")
        print(f"💡 Please add these keys to your policy.yaml file.")
        sys.exit(1)
    
    # Apply filters
    filtered_df = apply_guide_selection_filters(combined_df, CONFIG)
    
    if len(filtered_df) == 0:
        print("❌ No sequences passed the selection criteria")
        return
    
    # Select top guides
    top_guides = select_top_guides(filtered_df, CONFIG)
    
    # Save results
    output_file = CONFIG.get("OUTPUTS_TOP_GUIDES")
    top_guides.to_csv(output_file, index=False)
    print(f"✅ Saved {len(top_guides)} top guides to {output_file}")
    
    # Show summary
    print(f"\n📈 Selection Summary:")
    print(f"   Total sequences processed: {len(combined_df)}")
    print(f"   Sequences passing filters: {len(filtered_df)}")
    print(f"   Top guides selected: {len(top_guides)}")
    
    if len(top_guides) > 0:
        print(f"   Best combined score: {top_guides['on_plus_off'].max():.1f}")
        print(f"   Average combined score: {top_guides['on_plus_off'].mean():.1f}")
    
    print(f"\n💡 Open {output_file} in Excel and sort by 'on_plus_off' column to see best targets first!")

if __name__ == "__main__":
    main()
