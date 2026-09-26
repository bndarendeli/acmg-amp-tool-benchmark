#!/usr/bin/env python3
"""
FOXL2 Total Criterion Jaccard Analysis

Adds Total Criterion Jaccard to the existing FOXL2 capability-adjusted analysis.
Uses EXACTLY the same population and evaluability rules as the existing FOXL2
Mean Criterion Jaccard analysis.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Add path for criteria parser
sys.path.insert(0, str(Path(__file__).parent.parent / 'fig2'))
from utils_criteria_parser import parse_criteria_string

# Configuration
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'outputs' / 'foxl2_supplementary'

# Input files
CAPABILITY_MATRIX = BASE_DIR / 'outputs' / 'tool_criterion_capability_matrix.csv'
MASTER_TABLE = BASE_DIR / 'outputs' / 'fig4_master_variant_tool.csv'
EXISTING_SUMMARY = OUTPUT_DIR / 'foxl2_capability_adjusted_summary.csv'

# Output file
OUTPUT_FILE = OUTPUT_DIR / 'foxl2_total_criterion_jaccard.csv'


def load_capability_matrix():
    """Load operational criterion sets per tool"""
    print("\nLoading capability matrix...")
    df = pd.read_csv(CAPABILITY_MATRIX)
    
    operational_sets = {}
    for tool in df['tool'].unique():
        df_tool = df[df['tool'] == tool]
        operational = set(df_tool[df_tool['operational_status'] == 'supported']['criterion'])
        operational_sets[tool] = operational
        print(f"  {tool}: {len(operational)} operational criteria")
    
    return operational_sets


def calculate_raw_total_criterion_jaccard(df_tool):
    """
    Calculate raw Total Criterion Jaccard for a tool.
    
    Treats each (variant_id, criterion) assignment as a distinct element.
    """
    reference_assignments = []
    tool_assignments = []
    
    for _, row in df_tool.iterrows():
        variant_id = row['variant_id']
        
        # Parse reference criteria
        ref_crit = row['reference_criteria']
        if pd.notna(ref_crit):
            ref_set = parse_criteria_string(ref_crit)
            for crit in ref_set:
                reference_assignments.append((variant_id, crit))
        
        # Parse tool criteria
        tool_crit = row['tool_criteria']
        if pd.notna(tool_crit):
            tool_set = parse_criteria_string(tool_crit)
            for crit in tool_set:
                tool_assignments.append((variant_id, crit))
    
    # Convert to sets for Jaccard calculation
    ref_set = set(reference_assignments)
    tool_set = set(tool_assignments)
    
    # Calculate counts
    n_ref = len(ref_set)
    n_tool = len(tool_set)
    n_shared = len(ref_set.intersection(tool_set))
    n_fn = n_ref - n_shared
    n_fp = n_tool - n_shared
    
    # Calculate Jaccard
    union_size = n_ref + n_tool - n_shared
    if union_size == 0:
        total_jaccard = 1.0
    else:
        total_jaccard = n_shared / union_size
    
    return {
        'raw_reference_assignments': n_ref,
        'raw_tool_assignments': n_tool,
        'raw_shared_assignments': n_shared,
        'raw_FN': n_fn,
        'raw_FP': n_fp,
        'raw_total_criterion_jaccard': total_jaccard
    }


def calculate_capability_adjusted_total_criterion_jaccard(df_tool, operational_criteria):
    """
    Calculate capability-adjusted Total Criterion Jaccard for a tool.
    
    Projects criteria to operational scope, then pools (variant_id, criterion) assignments.
    Follows EXACTLY the same evaluability policy as existing FOXL2 analysis.
    """
    adjusted_reference_assignments = []
    adjusted_tool_assignments = []
    n_excluded = 0
    
    for _, row in df_tool.iterrows():
        variant_id = row['variant_id']
        
        # Parse criteria
        ref_crit = row['reference_criteria']
        tool_crit = row['tool_criteria']
        
        if pd.isna(ref_crit):
            ref_set = set()
        else:
            ref_set = parse_criteria_string(ref_crit)
        
        if pd.isna(tool_crit):
            tool_set = set()
        else:
            tool_set = parse_criteria_string(tool_crit)
        
        # Project to operational scope
        ref_supported = ref_set.intersection(operational_criteria)
        tool_supported = tool_set.intersection(operational_criteria)
        
        # FOXL2 evaluability policy: exclude if both projected sets are empty
        if len(ref_supported) == 0 and len(tool_supported) == 0:
            n_excluded += 1
            continue
        
        # Add projected assignments
        for crit in ref_supported:
            adjusted_reference_assignments.append((variant_id, crit))
        
        for crit in tool_supported:
            adjusted_tool_assignments.append((variant_id, crit))
    
    # Convert to sets for Jaccard calculation
    ref_set = set(adjusted_reference_assignments)
    tool_set = set(adjusted_tool_assignments)
    
    # Calculate counts
    n_ref = len(ref_set)
    n_tool = len(tool_set)
    n_shared = len(ref_set.intersection(tool_set))
    n_fn = n_ref - n_shared
    n_fp = n_tool - n_shared
    
    # Calculate Jaccard
    union_size = n_ref + n_tool - n_shared
    if union_size == 0:
        total_jaccard = 1.0
    else:
        total_jaccard = n_shared / union_size
    
    return {
        'adjusted_reference_assignments': n_ref,
        'adjusted_tool_assignments': n_tool,
        'adjusted_shared_assignments': n_shared,
        'adjusted_FN': n_fn,
        'adjusted_FP': n_fp,
        'adjusted_total_criterion_jaccard': total_jaccard,
        'n_excluded_both_empty': n_excluded
    }


def validate_metrics(raw_metrics, adj_metrics):
    """Validate mathematical consistency of metrics"""
    errors = []
    
    # Raw validation
    if raw_metrics['raw_shared_assignments'] > raw_metrics['raw_reference_assignments']:
        errors.append("raw_shared > raw_reference")
    
    if raw_metrics['raw_shared_assignments'] > raw_metrics['raw_tool_assignments']:
        errors.append("raw_shared > raw_tool")
    
    expected_ref = raw_metrics['raw_shared_assignments'] + raw_metrics['raw_FN']
    if expected_ref != raw_metrics['raw_reference_assignments']:
        errors.append(f"raw_reference != raw_shared + raw_FN ({raw_metrics['raw_reference_assignments']} != {expected_ref})")
    
    expected_tool = raw_metrics['raw_shared_assignments'] + raw_metrics['raw_FP']
    if expected_tool != raw_metrics['raw_tool_assignments']:
        errors.append(f"raw_tool != raw_shared + raw_FP ({raw_metrics['raw_tool_assignments']} != {expected_tool})")
    
    # Adjusted validation
    if adj_metrics['adjusted_shared_assignments'] > adj_metrics['adjusted_reference_assignments']:
        errors.append("adjusted_shared > adjusted_reference")
    
    if adj_metrics['adjusted_shared_assignments'] > adj_metrics['adjusted_tool_assignments']:
        errors.append("adjusted_shared > adjusted_tool")
    
    expected_ref = adj_metrics['adjusted_shared_assignments'] + adj_metrics['adjusted_FN']
    if expected_ref != adj_metrics['adjusted_reference_assignments']:
        errors.append(f"adjusted_reference != adjusted_shared + adjusted_FN ({adj_metrics['adjusted_reference_assignments']} != {expected_ref})")
    
    expected_tool = adj_metrics['adjusted_shared_assignments'] + adj_metrics['adjusted_FP']
    if expected_tool != adj_metrics['adjusted_tool_assignments']:
        errors.append(f"adjusted_tool != adjusted_shared + adjusted_FP ({adj_metrics['adjusted_tool_assignments']} != {expected_tool})")
    
    return errors


def main():
    """Main execution"""
    print("="*80)
    print("FOXL2 TOTAL CRITERION JACCARD ANALYSIS")
    print("="*80)
    print()
    print("Adding Total Criterion Jaccard to existing FOXL2 analysis")
    print("Using EXACT same population and evaluability rules")
    print()
    
    # Load data
    print("Loading data...")
    df_master = pd.read_csv(MASTER_TABLE)
    operational_sets = load_capability_matrix()
    df_existing = pd.read_csv(EXISTING_SUMMARY)
    
    # Filter to FOXL2, classification-concordant, evidence-evaluable
    # EXACTLY matching the existing FOXL2 analysis
    df_foxl2 = df_master[df_master['dataset'] == 'foxl2'].copy()
    print(f"\n  FOXL2 predictions: {len(df_foxl2):,}")
    
    df_concordant = df_foxl2[df_foxl2['classification_match'] == True].copy()
    print(f"  Classification-concordant: {len(df_concordant):,}")
    
    # Evidence-evaluable: criterion_jaccard is not NA
    df_evaluable = df_concordant[df_concordant['criterion_jaccard'].notna()].copy()
    print(f"  Evidence-evaluable (criterion_jaccard notna): {len(df_evaluable):,}")
    
    print("\n" + "="*80)
    print("CALCULATING TOTAL CRITERION JACCARD")
    print("="*80)
    
    results = []
    validation_passed = True
    
    for tool in sorted(df_evaluable['tool'].unique()):
        print(f"\n{tool}:")
        
        df_tool = df_evaluable[df_evaluable['tool'] == tool]
        n_raw = len(df_tool)
        
        # Get existing metrics for validation
        existing = df_existing[df_existing['tool'] == tool].iloc[0]
        existing_n_evaluable = existing['n_evaluable']
        existing_n_adjusted = existing['n_adjusted_evaluable']
        existing_raw_mean_j = existing['raw_mean_jaccard']
        existing_adj_mean_j = existing['adjusted_mean_jaccard']
        
        # Calculate raw Total CJ
        raw_metrics = calculate_raw_total_criterion_jaccard(df_tool)
        
        # Calculate adjusted Total CJ
        operational_criteria = operational_sets.get(tool, set())
        adj_metrics = calculate_capability_adjusted_total_criterion_jaccard(df_tool, operational_criteria)
        
        # Validate population
        print(f"  Population validation:")
        print(f"    n_raw_evaluable: {n_raw} (expected: {existing_n_evaluable})")
        
        if n_raw != existing_n_evaluable:
            print(f"    ✗ MISMATCH in raw population!")
            validation_passed = False
        else:
            print(f"    ✓ Raw population matches")
        
        n_adjusted = n_raw - adj_metrics['n_excluded_both_empty']
        print(f"    n_adjusted_evaluable: {n_adjusted} (expected: {existing_n_adjusted})")
        
        if n_adjusted != existing_n_adjusted:
            print(f"    ✗ MISMATCH in adjusted population!")
            validation_passed = False
        else:
            print(f"    ✓ Adjusted population matches")
        
        # Validate Mean CJ recomputation
        recomputed_raw_mean_j = df_tool['criterion_jaccard'].mean()
        print(f"  Mean CJ validation:")
        print(f"    Recomputed raw Mean CJ: {recomputed_raw_mean_j:.4f}")
        print(f"    Existing raw Mean CJ: {existing_raw_mean_j:.4f}")
        
        if abs(recomputed_raw_mean_j - existing_raw_mean_j) > 1e-6:
            print(f"    ✗ MISMATCH in raw Mean CJ!")
            validation_passed = False
        else:
            print(f"    ✓ Raw Mean CJ matches")
        
        # Validate mathematical consistency
        print(f"  Mathematical validation:")
        errors = validate_metrics(raw_metrics, adj_metrics)
        if errors:
            print(f"    ✗ Validation errors:")
            for error in errors:
                print(f"      - {error}")
            validation_passed = False
        else:
            print(f"    ✓ All mathematical checks passed")
        
        # Calculate delta
        delta_total_cj = adj_metrics['adjusted_total_criterion_jaccard'] - raw_metrics['raw_total_criterion_jaccard']
        
        # Store results
        result = {
            'tool': tool,
            'n_raw_evaluable': n_raw,
            'raw_mean_jaccard': existing_raw_mean_j,
            'raw_reference_assignments': raw_metrics['raw_reference_assignments'],
            'raw_tool_assignments': raw_metrics['raw_tool_assignments'],
            'raw_shared_assignments': raw_metrics['raw_shared_assignments'],
            'raw_FN': raw_metrics['raw_FN'],
            'raw_FP': raw_metrics['raw_FP'],
            'raw_total_criterion_jaccard': raw_metrics['raw_total_criterion_jaccard'],
            'n_adjusted_evaluable': n_adjusted,
            'adjusted_mean_jaccard': existing_adj_mean_j,
            'adjusted_reference_assignments': adj_metrics['adjusted_reference_assignments'],
            'adjusted_tool_assignments': adj_metrics['adjusted_tool_assignments'],
            'adjusted_shared_assignments': adj_metrics['adjusted_shared_assignments'],
            'adjusted_FN': adj_metrics['adjusted_FN'],
            'adjusted_FP': adj_metrics['adjusted_FP'],
            'adjusted_total_criterion_jaccard': adj_metrics['adjusted_total_criterion_jaccard'],
            'delta_total_criterion_jaccard': delta_total_cj
        }
        
        results.append(result)
        
        print(f"  Raw Total CJ: {raw_metrics['raw_total_criterion_jaccard']:.4f}")
        print(f"  Adjusted Total CJ: {adj_metrics['adjusted_total_criterion_jaccard']:.4f}")
        print(f"  Delta: {delta_total_cj:+.4f}")
    
    # Create DataFrame
    df_results = pd.DataFrame(results)
    
    # Save results
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df_results.to_csv(OUTPUT_FILE, index=False)
    
    print("\n" + "="*80)
    print("VALIDATION SUMMARY TABLE")
    print("="*80)
    print()
    print(f"{'Tool':<20} {'N raw':>8} {'Raw Mean J':>12} {'Raw Total J':>12} {'N adj':>8} {'Adj Mean J':>12} {'Adj Total J':>12} {'Delta Total J':>14}")
    print("-" * 130)
    
    for _, row in df_results.iterrows():
        print(f"{row['tool']:<20} {row['n_raw_evaluable']:>8} {row['raw_mean_jaccard']:>12.4f} {row['raw_total_criterion_jaccard']:>12.4f} {row['n_adjusted_evaluable']:>8} {row['adjusted_mean_jaccard']:>12.4f} {row['adjusted_total_criterion_jaccard']:>12.4f} {row['delta_total_criterion_jaccard']:>14.4f}")
    
    print()
    print("="*80)
    print("VALIDATION RESULTS")
    print("="*80)
    print()
    
    if validation_passed:
        print("✓ All validation checks PASSED")
        print("  ✓ Raw Mean CJ population validation passed for all tools")
        print("  ✓ Adjusted-evaluable n validation passed for all tools")
        print("  ✓ All Total CJ mathematical checks passed")
    else:
        print("✗ VALIDATION FAILED - see errors above")
    
    print()
    print(f"✓ Saved: {OUTPUT_FILE}")
    print()
    print("="*80)
    print("✓ FOXL2 TOTAL CRITERION JACCARD ANALYSIS COMPLETE")
    print("="*80)
    print()


if __name__ == '__main__':
    main()
