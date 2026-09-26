#!/usr/bin/env python3
"""
Generate Jaccard Heatmaps for ClinGen Datasets Only

Creates heatmaps matching analysis_2 format:
- Only ClinGen datasets (clingen_28012026, clingen_cancerpredisposition, clingen_hearing_loss)
- Criterion on rows (vertical), Tool on columns (horizontal)
- Color scheme: RdYlBu_r (Blue-Yellow-Orange-Red)
- Sorted by average performance
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / 'fig2'))
from utils_criteria_parser import parse_criteria_string

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 9
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

TOOLS = [
    'InterVar_2018', 'InterVar_2025', 'BIAS', 'CharGer_Local', 'CharGer_Online',
    'DiabloACMG', 'Exomiser', 'Genebe', 'TAPES', 'Franklin', 'AutoGVP', 'VIP-HL',
    'CancerSIGVAR', 'CPSR'
]

ALL_CRITERIA = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]


def is_criteria_available(criteria_str):
    """
    Determine if criteria data is available (not missing).
    
    Missing data indicators: 
    - None, NaN, empty string
    - 'Unknown', 'Not_Provided', 'NA', 'nan'
    - '.' (VCF-style missing value placeholder)
    
    CRITICAL: '.' in ClinGen data represents MISSING/UNAVAILABLE criterion
    information, NOT a genuinely observed empty criterion set.
    """
    if pd.isna(criteria_str):
        return False
    
    criteria_str = str(criteria_str).strip()
    if criteria_str in ['', '.', 'Unknown', 'Not_Provided', 'None', 'nan', 'NA']:
        return False
    
    return True


def calculate_jaccard_per_tool_criterion(df, dataset_name):
    """
    Calculate Jaccard similarity per tool per criterion.
    
    For each tool and criterion:
    - Find all variants where GT has the criterion
    - Find all variants where Tool has the criterion
    - Calculate Jaccard = |GT ∩ Tool| / |GT ∪ Tool|
    
    IMPORTANT: Excludes variants with missing/unavailable reference criteria
    (Ground_Truth_ACMG = '.') to avoid treating missing data as empty sets.
    
    Returns DataFrame with columns: Tool, Criterion, Dataset, Jaccard_Similarity, Intersection, Union
    """
    results = []
    
    for tool in TOOLS:
        criteria_col = f'{tool}_ACMG_Criteria'
        if criteria_col not in df.columns:
            continue
        
        for criterion in ALL_CRITERIA:
            gt_has_criterion = []
            tool_has_criterion = []
            
            for idx, row in df.iterrows():
                # CRITICAL FIX: Skip variants with missing reference criteria
                gt_acmg = row.get('Ground_Truth_ACMG', '')
                if not is_criteria_available(gt_acmg):
                    continue
                
                gt_criteria = parse_criteria_string(gt_acmg)
                tool_criteria = parse_criteria_string(row.get(criteria_col, ''))
                
                gt_has = criterion in gt_criteria
                tool_has = criterion in tool_criteria
                
                gt_has_criterion.append(gt_has)
                tool_has_criterion.append(tool_has)
            
            # Create sets of variant indices
            gt_set = set([i for i, x in enumerate(gt_has_criterion) if x])
            tool_set = set([i for i, x in enumerate(tool_has_criterion) if x])
            
            intersection = len(gt_set.intersection(tool_set))
            union = len(gt_set.union(tool_set))
            jaccard = intersection / union if union > 0 else 0.0
            
            results.append({
                'Tool': tool,
                'Criterion': criterion,
                'Dataset': dataset_name,
                'Jaccard_Similarity': jaccard,
                'Intersection': intersection,
                'Union': union
            })
    
    return pd.DataFrame(results)


def create_jaccard_heatmap(jaccard_df, dataset_name, output_dir):
    """
    Create Jaccard similarity heatmap matching analysis_2 format.
    
    Format:
    - Rows (index): ACMG Criteria (vertical axis)
    - Columns: Tools (horizontal axis)
    - Color: RdYlBu_r (Blue=high, Yellow=mid, Orange-Red=low)
    - Sorted by average performance
    """
    
    # Pivot: Criterion as rows (index), Tool as columns
    pivot_data = jaccard_df.pivot(index='Criterion', columns='Tool', values='Jaccard_Similarity')
    
    # Sort criteria by average performance (descending)
    pivot_data['avg'] = pivot_data.mean(axis=1)
    pivot_data = pivot_data.sort_values('avg', ascending=False)
    pivot_data = pivot_data.drop('avg', axis=1)
    
    # Create figure - height based on number of criteria
    fig, ax = plt.subplots(figsize=(14, max(8, len(pivot_data) * 0.3)))
    
    # Create heatmap with RdYlBu_r colormap (Blue-Yellow-Orange-Red)
    # Blue = high (1.0), Yellow = mid (0.5), Red = low (0.0)
    sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='RdYlBu_r', 
                vmin=0, vmax=1, cbar_kws={'label': 'Jaccard Similarity'}, 
                ax=ax, linewidths=0.5, linecolor='white')
    
    ax.set_title(f'Jaccard Similarity Heatmap - {dataset_name}', 
                 fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel('Tool', fontsize=10, fontweight='bold')
    ax.set_ylabel('ACMG Criterion', fontsize=10, fontweight='bold')
    
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(rotation=0, fontsize=9)
    
    plt.tight_layout()
    
    # Save with clean filename
    safe_name = dataset_name.replace(' ', '_').replace('/', '_').replace('+', '')
    plt.savefig(output_dir / f'jaccard_heatmap_{safe_name}.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / f'jaccard_heatmap_{safe_name}.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Created: jaccard_heatmap_{safe_name}.png/pdf")
    return safe_name


def load_dataset_with_tools(dataset_name, master_merged, gt_dir):
    """Load dataset ground truth and merge with tool results
    
    This matches the analysis_2 approach:
    1. Load ground truth Excel file
    2. Create Variant_Key
    3. Merge with master merged results
    
    Handles both ClinGen format (Chr, Pos, Ref, Alt) and FOXL2 format (hg38)
    """
    gt_file = gt_dir / f'{dataset_name}_ground_truth.xlsx'
    
    if not gt_file.exists():
        print(f"  ⚠ Ground truth file not found: {gt_file}")
        return None
    
    gt_df = pd.read_excel(gt_file)
    print(f"  Loaded {len(gt_df)} variants from ground truth")
    
    # Create variant key based on dataset format
    if 'hg38' in gt_df.columns:
        # FOXL2 format: hg38 column already has Variant_Key format
        gt_df['Variant_Key'] = gt_df['hg38'].astype(str)
        print(f"  Using hg38 column for Variant_Key (FOXL2 format)")
    elif 'Chr' in gt_df.columns:
        # ClinGen format: create from Chr-Pos-Ref-Alt
        gt_df['Variant_Key'] = (gt_df['Chr'].astype(str) + '-' + 
                                gt_df['Pos'].astype(str) + '-' + 
                                gt_df['Ref'].astype(str) + '-' + 
                                gt_df['Alt'].astype(str))
        print(f"  Created Variant_Key from Chr-Pos-Ref-Alt (ClinGen format)")
    else:
        print(f"  ⚠ Cannot create Variant_Key - missing required columns")
        return None
    
    # Merge with tool results
    merged = gt_df.merge(
        master_merged,
        on='Variant_Key',
        how='left',
        suffixes=('', '_master')
    )
    
    print(f"  Merged with tool results: {len(merged)} variants")
    
    return merged


def analyze_clingen_dataset(dataset_name, master_merged, gt_dir, results_dir, figures_dir):
    """Analyze a single ClinGen dataset"""
    
    print(f"\n{'='*80}")
    print(f"Analyzing: {dataset_name}")
    print(f"{'='*80}")
    
    # Load ground truth and merge with master results
    df = load_dataset_with_tools(dataset_name, master_merged, gt_dir)
    
    if df is None:
        return None
    
    print(f"  Calculating Jaccard similarity per criterion...")
    jaccard_df = calculate_jaccard_per_tool_criterion(df, dataset_name)
    
    print(f"  Creating Jaccard heatmap...")
    safe_name = create_jaccard_heatmap(jaccard_df, dataset_name, figures_dir)
    
    # Save data
    jaccard_df.to_csv(results_dir / f'jaccard_per_criterion_{safe_name}.csv', index=False)
    
    print(f"  ✓ Analysis complete for {dataset_name}")
    
    return jaccard_df, safe_name


def clean_old_figures(figures_dir, keep_files):
    """Delete old/incorrect heatmap figures"""
    print(f"\n{'='*80}")
    print("Cleaning old figures...")
    print(f"{'='*80}")
    
    # Find all jaccard heatmap files
    all_heatmaps = list(figures_dir.glob('jaccard_heatmap_*.png')) + \
                   list(figures_dir.glob('jaccard_heatmap_*.pdf'))
    
    deleted_count = 0
    for file in all_heatmaps:
        # Check if this file should be kept
        should_keep = any(keep_name in file.name for keep_name in keep_files)
        
        if not should_keep:
            print(f"  🗑️  Deleting: {file.name}")
            file.unlink()
            deleted_count += 1
    
    if deleted_count > 0:
        print(f"\n  ✓ Deleted {deleted_count} old heatmap files")
    else:
        print(f"\n  ✓ No old files to delete")


def main():
    """Main execution function"""
    
    base_dir = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last')
    analysis_2_dir = base_dir.parent / 'analysis_2'
    
    results_dir = base_dir / 'results'
    figures_dir = base_dir / 'figures'
    gt_dir = analysis_2_dir / 'ground_truth'
    
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("JACCARD HEATMAP GENERATION - CLINGEN DATASETS ONLY")
    print("="*80)
    print("\nColor Scheme: RdYlBu_r (Blue=High, Yellow=Mid, Red=Low)")
    print("Format: Criterion (rows) × Tool (columns)")
    print("\nMethod: Load ground truth + merge with master results (same as analysis_2)")
    
    # Load master merged results
    print(f"\n{'='*80}")
    print("Loading master merged results...")
    print(f"{'='*80}")
    master_file = analysis_2_dir / 'results' / 'all_variants' / 'merged_results.tsv'
    
    if not master_file.exists():
        print(f"\n✗ Master merged results not found: {master_file}")
        return
    
    master_merged = pd.read_csv(master_file, sep='\t')
    print(f"  Loaded {len(master_merged):,} variants from master file")
    
    # ClinGen datasets only (FOXL2 already correct)
    clingen_datasets = [
        'clingen_28012026',
    ]
    
    all_jaccard_results = []
    generated_files = []
    
    for dataset_name in clingen_datasets:
        result = analyze_clingen_dataset(dataset_name, master_merged, gt_dir, results_dir, figures_dir)
        
        if result is not None:
            jaccard_df, safe_name = result
            all_jaccard_results.append(jaccard_df)
            generated_files.append(safe_name)
    
    if all_jaccard_results:
        # Save combined results
        print(f"\n{'='*80}")
        print("Saving Combined Results")
        print(f"{'='*80}")
        
        combined_jaccard = pd.concat(all_jaccard_results, ignore_index=True)
        combined_jaccard.to_csv(results_dir / 'jaccard_clingen_datasets_combined.csv', index=False)
        print(f"  ✓ Saved: jaccard_clingen_datasets_combined.csv")
        
        # Clean old figures
        clean_old_figures(figures_dir, generated_files)
        
        print(f"\n{'='*80}")
        print("ANALYSIS COMPLETE")
        print(f"{'='*80}")
        print(f"\nGenerated {len(generated_files)} ClinGen heatmaps:")
        for name in generated_files:
            print(f"  - jaccard_heatmap_{name}.png/pdf")
        print(f"\nResults saved to: {results_dir}")
        print(f"Figures saved to: {figures_dir}")
        print("\n✓ All ClinGen Jaccard heatmaps generated successfully!")


if __name__ == "__main__":
    main()
