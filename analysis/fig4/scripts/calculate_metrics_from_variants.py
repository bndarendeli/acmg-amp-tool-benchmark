#!/usr/bin/env python3
"""
Calculate Figure 4 metrics directly from variant-level data.

Uses merged_results.tsv for ClinGen, HGMD–ClinVar HL, and HGMD–ClinVar Cancer.
Uses foxl2_merged_results.tsv for FOXL2.

Performs variant-level classification comparison with strict validation.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

print("=" * 80)
print("CALCULATING METRICS FROM VARIANT-LEVEL DATA")
print("=" * 80)

# Paths
base_dir = Path(__file__).parent.parent
data_dir = base_dir / 'data' / 'corrected'
data_dir.mkdir(parents=True, exist_ok=True)

source_dir = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last/data')
merged_results_file = source_dir / 'merged_results.tsv'
foxl2_results_file = source_dir / 'foxl2_merged_results.tsv'

# Cohort definitions
COHORTS = {
    'ClinGen': {
        'file': merged_results_file,
        'dataset_filter': 'clingen_28012026',
        'expected_n': 11409,
        'output_file': 'ClinGen_corrected_metrics.csv'
    },
    'FOXL2': {
        'file': foxl2_results_file,
        'dataset_filter': None,  # Use all variants
        'expected_n': 288,
        'output_file': 'FOXL2_corrected_metrics.csv'
    },
    'HGMD–ClinVar HL': {
        'file': merged_results_file,
        'dataset_filter': 'HGMD_Clinvar_HL',
        'expected_n': 1948,
        'output_file': 'HGMD+ClinVar_HL_corrected_metrics.csv'
    },
    'HGMD–ClinVar Cancer': {
        'file': merged_results_file,
        'dataset_filter': 'HGMD_Clinvar_Cancer',
        'expected_n': 1788,
        'output_file': 'HGMD+ClinVar_Cancer_corrected_metrics.csv'
    }
}

# Three-class normalization mapping
THREE_CLASS_MAP = {
    # P/LP variants
    'Pathogenic': 'P/LP',
    'Likely_Pathogenic': 'P/LP',
    'P': 'P/LP',
    'LP': 'P/LP',
    # VUS variants
    'VUS': 'VUS',
    'Uncertain_Significance': 'VUS',
    # B/LB variants
    'Benign': 'B/LB',
    'Likely_Benign': 'B/LB',
    'B': 'B/LB',
    'LB': 'B/LB'
}

# Valid ACMG 5-tier classifications (for call rate)
VALID_ACMG_CLASSES = {
    'Pathogenic', 'Likely_Pathogenic', 
    'Uncertain_Significance', 'VUS',
    'Likely_Benign', 'Benign'
}


def normalize_classification(value):
    """
    Normalize classification to 3-class system.
    Returns 'P/LP', 'VUS', 'B/LB', or None for invalid/missing.
    """
    if pd.isna(value):
        return None
    
    value_str = str(value).strip()
    
    # Check for no-call indicators
    if value_str in ['', 'Unknown', 'NOT_AVAILABLE', 'NULL', 'None', 'nan']:
        return None
    
    # Map to 3-class
    return THREE_CLASS_MAP.get(value_str, None)


def is_valid_call(value):
    """
    Check if prediction is a valid ACMG 5-tier classification.
    """
    if pd.isna(value):
        return False
    
    value_str = str(value).strip()
    
    if value_str in ['', 'Unknown', 'NOT_AVAILABLE', 'NULL', 'None', 'nan']:
        return False
    
    return value_str in VALID_ACMG_CLASSES


def calculate_metrics_for_tool(ground_truth, predictions, cohort_size):
    """
    Calculate all metrics for a single tool.
    
    Args:
        ground_truth: Series of normalized ground truth (P/LP, VUS, B/LB)
        predictions: Series of normalized predictions (P/LP, VUS, B/LB, or None for no-call)
        cohort_size: Total cohort size for validation
        
    Returns:
        Dictionary of metrics
    """
    # Validate input
    assert len(ground_truth) == len(predictions), "Mismatched lengths"
    assert len(ground_truth) == cohort_size, f"Expected {cohort_size} variants, got {len(ground_truth)}"
    
    # Count ground truth class distribution
    gt_counts = ground_truth.value_counts()
    support_plp = gt_counts.get('P/LP', 0)
    support_vus = gt_counts.get('VUS', 0)
    support_blb = gt_counts.get('B/LB', 0)
    
    # Initialize confusion counts
    metrics = {
        'P/LP': {'TP': 0, 'FP': 0, 'FN': 0},
        'VUS': {'TP': 0, 'FP': 0, 'FN': 0},
        'B/LB': {'TP': 0, 'FP': 0, 'FN': 0}
    }
    
    valid_calls = 0
    
    # Calculate confusion matrix elements
    for gt, pred in zip(ground_truth, predictions):
        if pred is not None:
            valid_calls += 1
            
            # True positives
            if gt == pred:
                metrics[gt]['TP'] += 1
            else:
                # False positive for predicted class
                metrics[pred]['FP'] += 1
                # False negative for true class
                metrics[gt]['FN'] += 1
        else:
            # No-call = false negative for true class
            metrics[gt]['FN'] += 1
    
    # Calculate per-class metrics
    results = {}
    f1_scores = []
    weights = []
    
    for class_name in ['P/LP', 'VUS', 'B/LB']:
        tp = metrics[class_name]['TP']
        fp = metrics[class_name]['FP']
        fn = metrics[class_name]['FN']
        
        # Determine ground truth support for this class
        if class_name == 'P/LP':
            gt_support = support_plp
        elif class_name == 'VUS':
            gt_support = support_vus
        else:  # B/LB
            gt_support = support_blb
        
        # If ground truth class has zero support, metrics are undefined (NaN)
        if gt_support == 0:
            recall = np.nan
            precision = np.nan
            f1 = np.nan
        else:
            # Recall
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            
            # Precision
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            
            # F1
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
        
        results[f'{class_name} Recall'] = recall
        results[f'{class_name} Precision'] = precision
        results[f'{class_name} F1'] = f1
        
        # For weighted F1 - only include classes with support > 0
        if gt_support > 0:
            f1_scores.append(f1)
            weights.append(gt_support)
    
    # Weighted F1 (using true ground-truth class supports)
    if sum(weights) > 0:
        weighted_f1 = sum(f * w for f, w in zip(f1_scores, weights)) / sum(weights)
    else:
        weighted_f1 = 0.0
    
    results['Weighted F1 (3-class)'] = weighted_f1
    
    # Call Rate
    call_rate = valid_calls / cohort_size if cohort_size > 0 else 0.0
    results['Call Rate'] = call_rate
    
    return results


def process_cohort(cohort_name, cohort_info):
    """
    Process a single cohort and calculate metrics for all tools.
    """
    print(f"\n{'='*80}")
    print(f"Processing {cohort_name} (expected n={cohort_info['expected_n']})")
    print(f"{'='*80}")
    
    # Load data
    file_path = cohort_info['file']
    if not file_path.exists():
        print(f"  ✗ File not found: {file_path}")
        return None
    
    print(f"  Loading: {file_path}")
    df = pd.read_csv(file_path, sep='\t', low_memory=False)
    
    # Filter by dataset if needed
    if cohort_info['dataset_filter']:
        df = df[df['Dataset'] == cohort_info['dataset_filter']].copy()
        print(f"  Filtered to Dataset == '{cohort_info['dataset_filter']}'")
    
    # Validate variant count
    n_variants = len(df)
    expected_n = cohort_info['expected_n']
    
    print(f"  Variants: {n_variants}")
    
    if n_variants != expected_n:
        print(f"  ⚠ WARNING: Expected {expected_n} variants, found {n_variants}")
    
    # Check for duplicates
    if 'Variant_Key' in df.columns:
        n_unique = df['Variant_Key'].nunique()
        if n_unique != n_variants:
            print(f"  ⚠ WARNING: Found {n_variants - n_unique} duplicate Variant_Keys")
    
    # Normalize ground truth
    df['GT_3class'] = df['Ground_Truth_Classification'].apply(normalize_classification)
    
    # Assert all ground truth classifications are valid
    invalid_gt = df['GT_3class'].isna()
    n_invalid = invalid_gt.sum()
    
    if n_invalid > 0:
        invalid_values = df.loc[invalid_gt, 'Ground_Truth_Classification'].unique()
        raise ValueError(
            f"Found {n_invalid} variants with invalid ground truth classifications in {cohort_name}.\n"
            f"Invalid values: {invalid_values}\n"
            f"All ground truth classifications must be valid 3-class mappings."
        )
    
    df_valid = df.copy()
    n_valid = len(df_valid)
    print(f"  ✓ All {n_valid} variants have valid ground truth classifications")
    
    # Calculate ground truth distribution
    gt_dist = df_valid['GT_3class'].value_counts()
    support_plp = gt_dist.get('P/LP', 0)
    support_vus = gt_dist.get('VUS', 0)
    support_blb = gt_dist.get('B/LB', 0)
    
    print(f"\n  Ground Truth Distribution:")
    print(f"    P/LP: {support_plp} ({support_plp/n_valid*100:.1f}%)")
    print(f"    VUS:  {support_vus} ({support_vus/n_valid*100:.1f}%)")
    print(f"    B/LB: {support_blb} ({support_blb/n_valid*100:.1f}%)")
    print(f"    Total: {support_plp + support_vus + support_blb}")
    
    # Find tool columns
    tool_columns = [col for col in df_valid.columns if col.endswith('_Classification')]
    tool_columns = [col for col in tool_columns if col != 'Ground_Truth_Classification']
    
    print(f"\n  Found {len(tool_columns)} tool columns")
    
    # Calculate metrics for each tool
    results = []
    
    for tool_col in sorted(tool_columns):
        tool_name = tool_col.replace('_Classification', '')
        
        # Normalize predictions
        df_valid['Pred_3class'] = df_valid[tool_col].apply(normalize_classification)
        
        # Calculate metrics
        try:
            metrics = calculate_metrics_for_tool(
                df_valid['GT_3class'],
                df_valid['Pred_3class'],
                n_valid
            )
            
            metrics['Tool'] = tool_name
            results.append(metrics)
            
            # Print summary
            print(f"  ✓ {tool_name}: F1={metrics['Weighted F1 (3-class)']:.3f}, "
                  f"Call Rate={metrics['Call Rate']*100:.1f}%")
            
        except Exception as e:
            print(f"  ✗ {tool_name}: Error - {e}")
    
    # Create DataFrame
    if not results:
        print(f"  ✗ No metrics calculated")
        return None
    
    results_df = pd.DataFrame(results)
    
    # Reorder columns
    col_order = ['Tool', 'P/LP Recall', 'B/LB Recall', 'VUS Recall', 
                 'P/LP Precision', 'B/LB Precision', 'VUS Precision',
                 'P/LP F1', 'B/LB F1', 'VUS F1',
                 'Weighted F1 (3-class)', 'Call Rate']
    results_df = results_df[col_order]
    
    # Save
    output_file = data_dir / cohort_info['output_file']
    results_df.to_csv(output_file, index=False, float_format='%.4f')
    print(f"\n  ✓ Saved: {output_file}")
    
    return results_df


# Process all cohorts
print("\n" + "=" * 80)
print("PROCESSING ALL COHORTS")
print("=" * 80)

all_results = {}

for cohort_name, cohort_info in COHORTS.items():
    results_df = process_cohort(cohort_name, cohort_info)
    if results_df is not None:
        all_results[cohort_name] = results_df

# Validation summary
print("\n" + "=" * 80)
print("VALIDATION SUMMARY")
print("=" * 80)

for cohort_name, results_df in all_results.items():
    print(f"\n{cohort_name}:")
    print(f"  Tools: {len(results_df)}")
    print(f"  Call Rate range: {results_df['Call Rate'].min():.3f} - {results_df['Call Rate'].max():.3f}")
    
    # Check if all call rates are valid
    if (results_df['Call Rate'] < 0).any() or (results_df['Call Rate'] > 1).any():
        print(f"  ⚠ WARNING: Invalid call rates detected!")

# Print metric tables
print("\n" + "=" * 80)
print("METRIC TABLES")
print("=" * 80)

for cohort_name, results_df in all_results.items():
    print(f"\n{cohort_name}:")
    print("-" * 120)
    
    # Format for display
    display_df = results_df[['Tool', 'P/LP Recall', 'B/LB Recall', 'VUS Recall', 
                             'Call Rate', 'Weighted F1 (3-class)']].copy()
    
    # Sort by Weighted F1
    display_df = display_df.sort_values('Weighted F1 (3-class)', ascending=False)
    
    # Print table
    print(f"{'Tool':<20} {'P/LP Recall':>12} {'B/LB Recall':>12} {'VUS Recall':>12} "
          f"{'Call Rate':>12} {'Weighted F1':>12}")
    print("-" * 120)
    
    for _, row in display_df.iterrows():
        print(f"{row['Tool']:<20} {row['P/LP Recall']:>12.3f} {row['B/LB Recall']:>12.3f} "
              f"{row['VUS Recall']:>12.3f} {row['Call Rate']:>12.3f} "
              f"{row['Weighted F1 (3-class)']:>12.3f}")

print("\n" + "=" * 80)
print("✓ METRIC CALCULATION COMPLETE")
print("=" * 80)
print(f"\nCorrected metrics saved to: {data_dir}")
print(f"\nNext step: Run create_figure4_dotplot.py to generate the figure")
print("=" * 80)
