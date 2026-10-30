#!/usr/bin/env python3
"""
Calculate Total Criterion Jaccard for ClinGen Tools

For each tool, treats each (variant, criterion) assignment as a distinct item
and calculates the Jaccard similarity across all assignments.

Uses the evidence-evaluable population from fig6_master_variant_tool.csv.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Configuration
OUTPUT_DIR = Path(__file__).parent / 'outputs'
MASTER_FILE = OUTPUT_DIR / 'fig6_master_variant_tool.csv'
CAPABILITY_MATRIX = OUTPUT_DIR / 'tool_criterion_capability_matrix.csv'
OUTPUT_FILE = OUTPUT_DIR / 'clingen_total_criterion_jaccard.csv'

# Tool order
TOOLS = [
    'InterVar_2018', 'InterVar_2025', 'BIAS', 'CharGer_Local', 'CharGer_Online',
    'DiabloACMG', 'Exomiser', 'Genebe', 'TAPES', 'Franklin', 'AutoGVP', 'VIP-HL',
    'CancerSIGVAR', 'CPSR'
]


def load_capability_matrix():
    """
    Load the capability matrix and create operational criterion sets per tool.
    
    Returns:
        dict: {tool: set of operationally supported criteria}
    """
    df_cap = pd.read_csv(CAPABILITY_MATRIX)
    
    # Build operational criterion sets
    operational_sets = {}
    for tool in df_cap['tool'].unique():
        df_tool = df_cap[df_cap['tool'] == tool]
        supported = set(df_tool[df_tool['operational_status'] == 'supported']['criterion'])
        operational_sets[tool] = supported
    
    return operational_sets


def parse_criteria_set(criteria_str):
    """Parse comma-separated criteria string into a set"""
    if pd.isna(criteria_str) or criteria_str == '':
        return set()
    return set(str(criteria_str).split(','))


def calculate_total_criterion_jaccard(df_tool):
    """
    Calculate raw total criterion Jaccard for a single tool.
    
    Treats each (variant, criterion) assignment as a distinct item.
    """
    # Count total assignments
    reference_assignments = 0
    tool_assignments = 0
    shared_assignments = 0
    
    for _, row in df_tool.iterrows():
        ref_criteria = parse_criteria_set(row['reference_criteria'])
        tool_criteria = parse_criteria_set(row['tool_criteria'])
        
        # Each criterion assigned to this variant counts as one assignment
        reference_assignments += len(ref_criteria)
        tool_assignments += len(tool_criteria)
        
        # Shared: same criterion assigned to same variant by both
        shared = ref_criteria & tool_criteria
        shared_assignments += len(shared)
    
    # Calculate false negatives and false positives
    false_negative_assignments = reference_assignments - shared_assignments
    false_positive_assignments = tool_assignments - shared_assignments
    
    # Calculate Total Criterion Jaccard
    denominator = reference_assignments + tool_assignments - shared_assignments
    total_criterion_jaccard = shared_assignments / denominator if denominator > 0 else 0.0
    
    return {
        'reference_assignments': reference_assignments,
        'tool_assignments': tool_assignments,
        'shared_assignments': shared_assignments,
        'false_negative_assignments': false_negative_assignments,
        'false_positive_assignments': false_positive_assignments,
        'total_criterion_jaccard': total_criterion_jaccard
    }


def calculate_capability_adjusted_total_criterion_jaccard(df_tool, operational_criteria):
    """
    Calculate capability-adjusted total criterion Jaccard for a single tool.
    
    Projects both reference and tool criteria to operational criteria,
    then treats each (variant, criterion) assignment as a distinct item.
    
    Args:
        df_tool: DataFrame with variant-tool rows
        operational_criteria: Set of operationally supported criteria for this tool
    
    Returns:
        dict: Adjusted metrics including n_adjusted_evaluable
    """
    # Count total assignments after projection
    adjusted_reference_assignments = 0
    adjusted_tool_assignments = 0
    adjusted_shared_assignments = 0
    n_adjusted_evaluable = 0
    
    for _, row in df_tool.iterrows():
        ref_criteria = parse_criteria_set(row['reference_criteria'])
        tool_criteria = parse_criteria_set(row['tool_criteria'])
        
        # Project to operational criteria (same logic as existing capability adjustment)
        adjusted_ref = ref_criteria & operational_criteria
        adjusted_tool = tool_criteria & operational_criteria
        
        # Follow existing policy: if both projected sets are empty, skip this row
        # (same as adjusted_jaccard = None in perform_capability_adjusted_analysis.py)
        if len(adjusted_ref) == 0 and len(adjusted_tool) == 0:
            continue
        
        # This row has operationally evaluable criteria
        n_adjusted_evaluable += 1
        
        # Count assignments
        adjusted_reference_assignments += len(adjusted_ref)
        adjusted_tool_assignments += len(adjusted_tool)
        
        # Shared: same criterion assigned to same variant by both (after projection)
        shared = adjusted_ref & adjusted_tool
        adjusted_shared_assignments += len(shared)
    
    # Calculate false negatives and false positives
    adjusted_false_negative_assignments = adjusted_reference_assignments - adjusted_shared_assignments
    adjusted_false_positive_assignments = adjusted_tool_assignments - adjusted_shared_assignments
    
    # Calculate Adjusted Total Criterion Jaccard
    denominator = adjusted_reference_assignments + adjusted_tool_assignments - adjusted_shared_assignments
    adjusted_total_criterion_jaccard = adjusted_shared_assignments / denominator if denominator > 0 else 0.0
    
    return {
        'n_adjusted_evaluable': n_adjusted_evaluable,
        'adjusted_reference_assignments': adjusted_reference_assignments,
        'adjusted_tool_assignments': adjusted_tool_assignments,
        'adjusted_shared_assignments': adjusted_shared_assignments,
        'adjusted_false_negative_assignments': adjusted_false_negative_assignments,
        'adjusted_false_positive_assignments': adjusted_false_positive_assignments,
        'adjusted_total_criterion_jaccard': adjusted_total_criterion_jaccard
    }


def main():
    """Main execution"""
    print("="*80)
    print("CALCULATE TOTAL CRITERION JACCARD")
    print("="*80)
    print()
    
    # Load master table
    print(f"Loading master table: {MASTER_FILE}")
    df_master = pd.read_csv(MASTER_FILE)
    print(f"  ✓ Loaded {len(df_master):,} rows")
    
    # Filter to ClinGen evidence-evaluable variants
    # CRITICAL: Use EXACTLY the same rows as raw Mean Criterion Jaccard:
    # 1. ClinGen dataset
    # 2. Classification-concordant (classification_match == True)
    # 3. Evaluable criterion evidence (criterion_jaccard is not null)
    df_clingen = df_master[
        (df_master['dataset'] == 'clingen_28012026') &
        (df_master['classification_match'] == True) &
        (df_master['criterion_jaccard'].notna())
    ].copy()
    
    print(f"\nClinGen evidence-evaluable population (same as raw Mean Criterion Jaccard):")
    print(f"  Dataset: clingen_28012026")
    print(f"  Classification-concordant: classification_match == True")
    print(f"  Evaluable evidence: criterion_jaccard.notna()")
    print(f"  Variants × Tools: {len(df_clingen):,}")
    print(f"  Unique variants: {df_clingen['variant_id'].nunique():,}")
    print()
    
    # Load existing mean criterion jaccard for validation
    print("Loading existing Mean Criterion Jaccard for validation...")
    df_existing = pd.read_csv(OUTPUT_DIR / 'fig6_evidence_distributions.csv')
    df_existing = df_existing[
        (df_existing['dataset'] == 'clingen_28012026') &
        (df_existing['tool'] != 'POOLED')
    ].copy()
    print(f"  ✓ Loaded existing metrics for {len(df_existing)} tools")
    
    # Load existing capability-adjusted summary for validation
    print("\nLoading existing capability-adjusted summary for validation...")
    df_cap_summary = pd.read_csv(OUTPUT_DIR / 'fig6_capability_adjusted_summary.csv')
    print(f"  ✓ Loaded capability-adjusted metrics for {len(df_cap_summary)} tools")
    
    # Load capability matrix
    print("\nLoading capability matrix...")
    operational_sets = load_capability_matrix()
    print(f"  ✓ Loaded operational criteria for {len(operational_sets)} tools")
    print()
    
    # Calculate Total Criterion Jaccard for each tool
    print("Calculating Total Criterion Jaccard for each tool...")
    print()
    
    results = []
    
    for tool in TOOLS:
        df_tool = df_clingen[df_clingen['tool'] == tool].copy()
        n_evaluable = len(df_tool)
        
        if n_evaluable == 0:
            print(f"  ⚠ {tool}: No evaluable variants")
            continue
        
        # Recompute mean criterion jaccard from the same rows
        recomputed_mean_jaccard = df_tool['criterion_jaccard'].mean()
        
        # Get existing mean criterion jaccard
        existing_row = df_existing[df_existing['tool'] == tool]
        if len(existing_row) > 0:
            existing_mean_jaccard = existing_row.iloc[0]['mean_criterion_jaccard']
            mean_match = abs(recomputed_mean_jaccard - existing_mean_jaccard) < 0.0001
        else:
            existing_mean_jaccard = None
            mean_match = False
        
        # Get existing capability-adjusted n_adjusted_evaluable for validation
        cap_row = df_cap_summary[df_cap_summary['tool'] == tool]
        if len(cap_row) > 0:
            existing_n_adjusted = cap_row.iloc[0]['n_adjusted_evaluable']
        else:
            existing_n_adjusted = None
        
        # Calculate raw metrics
        raw_metrics = calculate_total_criterion_jaccard(df_tool)
        
        # Calculate capability-adjusted metrics
        operational_criteria = operational_sets.get(tool, set())
        adjusted_metrics = calculate_capability_adjusted_total_criterion_jaccard(df_tool, operational_criteria)
        
        # Check n_adjusted_evaluable match
        n_adjusted_match = (adjusted_metrics['n_adjusted_evaluable'] == existing_n_adjusted) if existing_n_adjusted is not None else False
        
        # Calculate delta
        delta_total_criterion_jaccard = adjusted_metrics['adjusted_total_criterion_jaccard'] - raw_metrics['total_criterion_jaccard']
        
        results.append({
            'tool': tool,
            'n_evaluable': n_evaluable,
            'recomputed_mean_jaccard': recomputed_mean_jaccard,
            'existing_mean_jaccard': existing_mean_jaccard,
            'mean_match': mean_match,
            'raw_reference_assignments': raw_metrics['reference_assignments'],
            'raw_tool_assignments': raw_metrics['tool_assignments'],
            'raw_shared_assignments': raw_metrics['shared_assignments'],
            'raw_FN': raw_metrics['false_negative_assignments'],
            'raw_FP': raw_metrics['false_positive_assignments'],
            'raw_total_criterion_jaccard': raw_metrics['total_criterion_jaccard'],
            'n_adjusted_evaluable': adjusted_metrics['n_adjusted_evaluable'],
            'existing_n_adjusted_evaluable': existing_n_adjusted,
            'n_adjusted_match': n_adjusted_match,
            'adjusted_reference_assignments': adjusted_metrics['adjusted_reference_assignments'],
            'adjusted_tool_assignments': adjusted_metrics['adjusted_tool_assignments'],
            'adjusted_shared_assignments': adjusted_metrics['adjusted_shared_assignments'],
            'adjusted_FN': adjusted_metrics['adjusted_false_negative_assignments'],
            'adjusted_FP': adjusted_metrics['adjusted_false_positive_assignments'],
            'adjusted_total_criterion_jaccard': adjusted_metrics['adjusted_total_criterion_jaccard'],
            'delta_total_criterion_jaccard': delta_total_criterion_jaccard
        })
        
        match_str = "✓" if mean_match else "✗"
        n_adj_str = "✓" if n_adjusted_match else "✗"
        print(f"  {match_str} {tool}: n={n_evaluable:,}, n_adj={adjusted_metrics['n_adjusted_evaluable']:,}{n_adj_str}, "
              f"Raw Total J={raw_metrics['total_criterion_jaccard']:.4f}, "
              f"Adj Total J={adjusted_metrics['adjusted_total_criterion_jaccard']:.4f}, "
              f"Δ={delta_total_criterion_jaccard:+.4f}")
    
    # Create DataFrame
    df_results = pd.DataFrame(results)
    
    # Validation checks
    print()
    print("="*80)
    print("VALIDATION CHECKS")
    print("="*80)
    print()
    
    all_valid = True
    mean_validation_passed = True
    n_adjusted_validation_passed = True
    
    for _, row in df_results.iterrows():
        tool = row['tool']
        
        # Check 0: Mean Criterion Jaccard match
        if not row['mean_match']:
            print(f"  ✗ {tool}: Recomputed mean != existing mean")
            mean_validation_passed = False
            all_valid = False
        
        # Check 0b: n_adjusted_evaluable match
        if not row['n_adjusted_match']:
            print(f"  ✗ {tool}: n_adjusted_evaluable mismatch (got {row['n_adjusted_evaluable']}, expected {row['existing_n_adjusted_evaluable']})")
            n_adjusted_validation_passed = False
            all_valid = False
        
        # RAW METRICS VALIDATION
        # Check 1: raw shared <= raw reference
        if row['raw_shared_assignments'] > row['raw_reference_assignments']:
            print(f"  ✗ {tool}: raw shared > raw reference")
            all_valid = False
        
        # Check 2: raw shared <= raw tool
        if row['raw_shared_assignments'] > row['raw_tool_assignments']:
            print(f"  ✗ {tool}: raw shared > raw tool")
            all_valid = False
        
        # Check 3: raw reference = raw shared + raw FN
        expected_ref = row['raw_shared_assignments'] + row['raw_FN']
        if abs(row['raw_reference_assignments'] - expected_ref) > 0.001:
            print(f"  ✗ {tool}: raw reference != raw shared + raw FN")
            all_valid = False
        
        # Check 4: raw tool = raw shared + raw FP
        expected_tool = row['raw_shared_assignments'] + row['raw_FP']
        if abs(row['raw_tool_assignments'] - expected_tool) > 0.001:
            print(f"  ✗ {tool}: raw tool != raw shared + raw FP")
            all_valid = False
        
        # Check 5: raw Jaccard calculation
        denominator = row['raw_reference_assignments'] + row['raw_tool_assignments'] - row['raw_shared_assignments']
        expected_jaccard = row['raw_shared_assignments'] / denominator if denominator > 0 else 0.0
        if abs(row['raw_total_criterion_jaccard'] - expected_jaccard) > 0.0001:
            print(f"  ✗ {tool}: raw Jaccard calculation incorrect")
            all_valid = False
        
        # ADJUSTED METRICS VALIDATION
        # Check 6: adjusted shared <= adjusted reference
        if row['adjusted_shared_assignments'] > row['adjusted_reference_assignments']:
            print(f"  ✗ {tool}: adjusted shared > adjusted reference")
            all_valid = False
        
        # Check 7: adjusted shared <= adjusted tool
        if row['adjusted_shared_assignments'] > row['adjusted_tool_assignments']:
            print(f"  ✗ {tool}: adjusted shared > adjusted tool")
            all_valid = False
        
        # Check 8: adjusted reference = adjusted shared + adjusted FN
        expected_adj_ref = row['adjusted_shared_assignments'] + row['adjusted_FN']
        if abs(row['adjusted_reference_assignments'] - expected_adj_ref) > 0.001:
            print(f"  ✗ {tool}: adjusted reference != adjusted shared + adjusted FN")
            all_valid = False
        
        # Check 9: adjusted tool = adjusted shared + adjusted FP
        expected_adj_tool = row['adjusted_shared_assignments'] + row['adjusted_FP']
        if abs(row['adjusted_tool_assignments'] - expected_adj_tool) > 0.001:
            print(f"  ✗ {tool}: adjusted tool != adjusted shared + adjusted FP")
            all_valid = False
        
        # Check 10: adjusted Jaccard calculation
        adj_denominator = row['adjusted_reference_assignments'] + row['adjusted_tool_assignments'] - row['adjusted_shared_assignments']
        expected_adj_jaccard = row['adjusted_shared_assignments'] / adj_denominator if adj_denominator > 0 else 0.0
        if abs(row['adjusted_total_criterion_jaccard'] - expected_adj_jaccard) > 0.0001:
            print(f"  ✗ {tool}: adjusted Jaccard calculation incorrect")
            all_valid = False
    
    if all_valid:
        print("  ✓ All validation checks passed")
        if mean_validation_passed:
            print("  ✓ All recomputed means match existing means (population verified)")
        if n_adjusted_validation_passed:
            print("  ✓ All n_adjusted_evaluable match existing capability-adjusted analysis")
    else:
        print("  ⚠ Some validation checks failed")
    
    # Save results
    print()
    print("="*80)
    print("SAVING RESULTS")
    print("="*80)
    print()
    
    df_results.to_csv(OUTPUT_FILE, index=False, float_format='%.6f')
    print(f"  ✓ Saved: {OUTPUT_FILE}")
    
    # Print compact summary table
    print()
    print("="*80)
    print("TOTAL CRITERION JACCARD RESULTS (RAW vs ADJUSTED)")
    print("="*80)
    print()
    
    # Sort by Raw Total Criterion Jaccard descending
    df_display = df_results.sort_values('raw_total_criterion_jaccard', ascending=False)
    
    print(f"{'Tool':<20} {'Raw Total J':>12} {'Adj Total J':>12} {'Delta':>8} {'N raw':>8} {'N adj':>8}")
    print("-"*80)
    
    for _, row in df_display.iterrows():
        print(f"{row['tool']:<20} {row['raw_total_criterion_jaccard']:>12.6f} "
              f"{row['adjusted_total_criterion_jaccard']:>12.6f} "
              f"{row['delta_total_criterion_jaccard']:>+8.4f} "
              f"{row['n_evaluable']:>8.0f} {row['n_adjusted_evaluable']:>8.0f}")
    
    print()
    print("="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print()
    
    print(f"Tools: {len(df_results)}")
    
    print(f"\nRaw Total Criterion Jaccard:")
    print(f"  Mean: {df_results['raw_total_criterion_jaccard'].mean():.6f}")
    print(f"  Median: {df_results['raw_total_criterion_jaccard'].median():.6f}")
    print(f"  Range: [{df_results['raw_total_criterion_jaccard'].min():.6f}, {df_results['raw_total_criterion_jaccard'].max():.6f}]")
    
    print(f"\nAdjusted Total Criterion Jaccard:")
    print(f"  Mean: {df_results['adjusted_total_criterion_jaccard'].mean():.6f}")
    print(f"  Median: {df_results['adjusted_total_criterion_jaccard'].median():.6f}")
    print(f"  Range: [{df_results['adjusted_total_criterion_jaccard'].min():.6f}, {df_results['adjusted_total_criterion_jaccard'].max():.6f}]")
    
    print(f"\nDelta (Adjusted - Raw):")
    print(f"  Mean: {df_results['delta_total_criterion_jaccard'].mean():+.6f}")
    print(f"  Median: {df_results['delta_total_criterion_jaccard'].median():+.6f}")
    print(f"  Range: [{df_results['delta_total_criterion_jaccard'].min():+.6f}, {df_results['delta_total_criterion_jaccard'].max():+.6f}]")
    print()
    
    print("="*80)
    print("✓ TOTAL CRITERION JACCARD CALCULATION COMPLETE")
    print("="*80)


if __name__ == '__main__':
    main()
