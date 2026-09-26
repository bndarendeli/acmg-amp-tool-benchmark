#!/usr/bin/env python3
"""
Reconstruct Classification Performance Metrics for Figure 4

Calculates 3-class Accuracy and class-specific F1 scores using the complete
reference cohort as denominator (including variants without valid tool predictions).

Metrics:
- Accuracy (3-class): Correct predictions / Total reference cohort
- P/LP F1, VUS F1, B/LB F1: Class-specific F1 scores

Cohorts:
- ClinGen Master Cohort: n = 11,409
- FOXL2: n = 288
- HGMD–ClinVar Hearing Loss: n = 1,948
- HGMD–ClinVar Cancer Predisposition: n = 1,788
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

print("=" * 80)
print("RECONSTRUCTING CLASSIFICATION PERFORMANCE METRICS")
print("=" * 80)

# Paths
base_dir = Path(__file__).parent.parent
output_dir = base_dir / 'data' / 'reconstructed'
output_dir.mkdir(parents=True, exist_ok=True)

# Source data files
source_dir = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last/data')
merged_results_file = source_dir / 'merged_results.tsv'
foxl2_results_file = source_dir / 'foxl2_merged_results.tsv'

# Cohort definitions
COHORTS = {
    'ClinGen': {
        'file': merged_results_file,
        'dataset_filter': 'clingen_28012026',
        'expected_n': 11409,
        'output_file': 'ClinGen_reconstructed_metrics.csv'
    },
    'FOXL2': {
        'file': foxl2_results_file,
        'dataset_filter': None,
        'expected_n': 288,
        'output_file': 'FOXL2_reconstructed_metrics.csv'
    },
    'HGMD–ClinVar HL': {
        'file': merged_results_file,
        'dataset_filter': 'HGMD_Clinvar_HL',
        'expected_n': 1948,
        'output_file': 'HGMD_ClinVar_HL_reconstructed_metrics.csv'
    },
    'HGMD–ClinVar Cancer': {
        'file': merged_results_file,
        'dataset_filter': 'HGMD_Clinvar_Cancer',
        'expected_n': 1788,
        'output_file': 'HGMD_ClinVar_Cancer_reconstructed_metrics.csv'
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


def calculate_metrics(df, reference_col, tool_col, cohort_size):
    """
    Calculate 3-class Accuracy and class-specific F1 scores.
    
    Uses complete cohort as denominator:
    - Variants without valid tool predictions count as incorrect
    - For F1, they count as false negatives for their reference class
    
    Args:
        df: DataFrame with reference and tool classifications
        reference_col: Column name for reference classification
        tool_col: Column name for tool classification
        cohort_size: Total reference cohort size
    
    Returns:
        dict: Metrics including accuracy and class-specific F1 scores
    """
    # Normalize classifications
    df['ref_3class'] = df[reference_col].apply(normalize_classification)
    df['tool_3class'] = df[tool_col].apply(normalize_classification)
    
    # Count reference class support
    ref_counts = df['ref_3class'].value_counts()
    support_plp = ref_counts.get('P/LP', 0)
    support_vus = ref_counts.get('VUS', 0)
    support_blb = ref_counts.get('B/LB', 0)
    
    # Initialize confusion matrix
    metrics = {
        'P/LP': {'TP': 0, 'FP': 0, 'FN': 0},
        'VUS': {'TP': 0, 'FP': 0, 'FN': 0},
        'B/LB': {'TP': 0, 'FP': 0, 'FN': 0}
    }
    
    # Count correct predictions for accuracy
    n_correct = 0
    n_valid_predictions = 0
    
    # Process each variant
    for _, row in df.iterrows():
        ref_class = row['ref_3class']
        tool_class = row['tool_3class']
        
        # Skip if reference is invalid
        if ref_class is None:
            continue
        
        if tool_class is not None:
            # Valid tool prediction
            n_valid_predictions += 1
            
            if ref_class == tool_class:
                # Correct prediction
                n_correct += 1
                metrics[ref_class]['TP'] += 1
            else:
                # Incorrect prediction
                metrics[tool_class]['FP'] += 1
                metrics[ref_class]['FN'] += 1
        else:
            # No valid tool prediction = false negative for reference class
            metrics[ref_class]['FN'] += 1
    
    # Calculate 3-class Accuracy
    accuracy_3class = n_correct / cohort_size if cohort_size > 0 else 0.0
    
    # Calculate class-specific metrics
    results = {
        'cohort_n': cohort_size,
        'n_valid_predictions': n_valid_predictions,
        'n_correct_3class': n_correct,
        'accuracy_3class': accuracy_3class
    }
    
    for class_name in ['P/LP', 'VUS', 'B/LB']:
        tp = metrics[class_name]['TP']
        fp = metrics[class_name]['FP']
        fn = metrics[class_name]['FN']
        
        # Determine ground truth support
        if class_name == 'P/LP':
            gt_support = support_plp
        elif class_name == 'VUS':
            gt_support = support_vus
        else:  # B/LB
            gt_support = support_blb
        
        # If ground truth class has zero support, metrics are NA
        if gt_support == 0:
            precision = np.nan
            recall = np.nan
            f1 = np.nan
        else:
            # Precision
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            
            # Recall
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            
            # F1
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
        
        results[f'{class_name.replace("/", "_")}_precision'] = precision
        results[f'{class_name.replace("/", "_")}_recall'] = recall
        results[f'{class_name.replace("/", "_")}_f1'] = f1
        results[f'{class_name.replace("/", "_")}_support'] = gt_support
    
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
    if 'Variant_ID' in df.columns:
        n_duplicates = df['Variant_ID'].duplicated().sum()
        if n_duplicates > 0:
            print(f"  ⚠ WARNING: {n_duplicates} duplicate variant IDs found")
    
    # Identify reference column
    if 'Ground_Truth_Classification' in df.columns:
        reference_col = 'Ground_Truth_Classification'
    elif 'Ground_Truth' in df.columns:
        reference_col = 'Ground_Truth'
    elif 'Reference' in df.columns:
        reference_col = 'Reference'
    else:
        print(f"  ✗ ERROR: No reference column found")
        print(f"  Available columns: {df.columns.tolist()[:10]}")
        return None
    
    # Identify tool columns (those ending with _Classification)
    tool_columns = [col for col in df.columns if col.endswith('_Classification') 
                   and col != reference_col]
    
    print(f"  Tool columns: {len(tool_columns)}")
    
    # Check for unmapped classifications
    ref_values = df[reference_col].dropna().unique()
    unmapped_ref = [v for v in ref_values if normalize_classification(v) is None]
    if unmapped_ref:
        print(f"  ⚠ WARNING: Unmapped reference classifications: {unmapped_ref}")
    
    # Calculate metrics for each tool
    results = []
    
    for tool_col in tool_columns:
        print(f"    Processing {tool_col}...", end='')
        
        # Check for unmapped tool classifications
        tool_values = df[tool_col].dropna().unique()
        unmapped_tool = [v for v in tool_values if normalize_classification(v) is None 
                        and str(v) not in ['', 'Unknown', 'NOT_AVAILABLE', 'NULL', 'None', 'nan']]
        if unmapped_tool:
            print(f"\n      ⚠ WARNING: Unmapped tool classifications: {unmapped_tool}")
        
        metrics = calculate_metrics(df, reference_col, tool_col, expected_n)
        metrics['dataset'] = cohort_name
        # Remove "_Classification" suffix from tool name
        metrics['tool'] = tool_col.replace('_Classification', '')
        results.append(metrics)
        print(f" ✓")
    
    # Create DataFrame
    df_results = pd.DataFrame(results)
    
    # Reorder columns
    column_order = [
        'dataset', 'tool', 'cohort_n', 'n_valid_predictions', 'n_correct_3class', 
        'accuracy_3class',
        'P_LP_precision', 'P_LP_recall', 'P_LP_f1', 'P_LP_support',
        'VUS_precision', 'VUS_recall', 'VUS_f1', 'VUS_support',
        'B_LB_precision', 'B_LB_recall', 'B_LB_f1', 'B_LB_support'
    ]
    df_results = df_results[column_order]
    
    # Save
    output_file = output_dir / cohort_info['output_file']
    df_results.to_csv(output_file, index=False)
    print(f"  ✓ Saved: {output_file}")
    
    # Validation
    print(f"\n  Validation:")
    print(f"    Tool instances: {len(df_results)}")
    print(f"    Cohort denominator check: All rows have cohort_n = {expected_n}: {(df_results['cohort_n'] == expected_n).all()}")
    print(f"    Accuracy uses cohort_n: {(df_results['n_correct_3class'] <= df_results['cohort_n']).all()}")
    
    return df_results


def main():
    """Main execution"""
    
    all_results = {}
    
    for cohort_name, cohort_info in COHORTS.items():
        df_results = process_cohort(cohort_name, cohort_info)
        if df_results is not None:
            all_results[cohort_name] = df_results
    
    # Print ClinGen summary table
    if 'ClinGen' in all_results:
        print(f"\n{'='*80}")
        print("CLINGEN MASTER COHORT SUMMARY")
        print(f"{'='*80}\n")
        
        df_clingen = all_results['ClinGen']
        
        # Format for display
        display_cols = ['tool', 'cohort_n', 'n_valid_predictions', 'n_correct_3class', 
                       'accuracy_3class', 'P_LP_f1', 'VUS_f1', 'B_LB_f1']
        df_display = df_clingen[display_cols].copy()
        
        # Format numeric columns
        df_display['accuracy_3class'] = df_display['accuracy_3class'].apply(lambda x: f'{x:.4f}')
        df_display['P_LP_f1'] = df_display['P_LP_f1'].apply(lambda x: f'{x:.4f}' if not pd.isna(x) else 'NA')
        df_display['VUS_f1'] = df_display['VUS_f1'].apply(lambda x: f'{x:.4f}' if not pd.isna(x) else 'NA')
        df_display['B_LB_f1'] = df_display['B_LB_f1'].apply(lambda x: f'{x:.4f}' if not pd.isna(x) else 'NA')
        
        print(df_display.to_string(index=False))
    
    # Summary
    print(f"\n{'='*80}")
    print("RECONSTRUCTION COMPLETE")
    print(f"{'='*80}\n")
    
    print("Source Files Used:")
    print(f"  - {merged_results_file}")
    print(f"  - {foxl2_results_file}")
    
    print(f"\nFiles Created:")
    for cohort_name, cohort_info in COHORTS.items():
        output_file = output_dir / cohort_info['output_file']
        if output_file.exists():
            print(f"  ✓ {output_file}")
    
    print(f"\n{'='*80}")
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
