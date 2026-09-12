#!/usr/bin/env python3
"""
Generate FOXL2 Criterion-Level Missing Decomposition

Creates criterion-level decomposition from variant-level capability-adjusted data.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Configuration
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'outputs' / 'foxl2_supplementary'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Input file
VARIANT_LEVEL_FILE = OUTPUT_DIR / 'foxl2_capability_adjusted_variant_tool.csv'

# Output file
DECOMP_FILE = OUTPUT_DIR / 'foxl2_criterion_missing_decomposition.csv'

# Canonical 28 criteria
CANONICAL_CRITERIA = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]


def parse_criteria_list(criteria_str):
    """Parse comma-separated criteria string into set"""
    if pd.isna(criteria_str) or criteria_str == '':
        return set()
    return set(str(criteria_str).split(','))


def generate_criterion_decomposition():
    """Generate criterion-level missing decomposition from variant-level data"""
    print("="*80)
    print("GENERATING FOXL2 CRITERION-LEVEL DECOMPOSITION")
    print("="*80)
    print()
    
    # Load variant-level data
    print(f"Loading: {VARIANT_LEVEL_FILE.name}")
    df = pd.read_csv(VARIANT_LEVEL_FILE)
    print(f"  Total variant-tool pairs: {len(df):,}")
    
    # Initialize counters for each criterion
    criterion_stats = {crit: {
        'n_outside_scope': 0,
        'n_supported_but_missed': 0
    } for crit in CANONICAL_CRITERIA}
    
    # Process each variant-tool pair
    print("\nProcessing variant-tool pairs...")
    for idx, row in df.iterrows():
        # Parse missing criteria
        outside_scope = parse_criteria_list(row['outside_scope_missing_criteria'])
        supported_missed = parse_criteria_list(row['supported_but_missed_criteria'])
        
        # Count occurrences
        for crit in outside_scope:
            if crit in criterion_stats:
                criterion_stats[crit]['n_outside_scope'] += 1
        
        for crit in supported_missed:
            if crit in criterion_stats:
                criterion_stats[crit]['n_supported_but_missed'] += 1
    
    # Build results DataFrame
    results = []
    for crit in CANONICAL_CRITERIA:
        n_outside = criterion_stats[crit]['n_outside_scope']
        n_supported = criterion_stats[crit]['n_supported_but_missed']
        n_total = n_outside + n_supported
        
        # Calculate fractions
        if n_total > 0:
            outside_frac = n_outside / n_total
            supported_frac = n_supported / n_total
        else:
            outside_frac = 0.0
            supported_frac = 0.0
        
        results.append({
            'criterion': crit,
            'n_total_missing': n_total,
            'n_outside_operational_scope': n_outside,
            'n_operationally_evaluable_but_not_assigned': n_supported,
            'outside_scope_fraction_of_missing': outside_frac,
            'operationally_evaluable_but_not_assigned_fraction': supported_frac
        })
    
    df_results = pd.DataFrame(results)
    
    # Validation
    print("\nValidation:")
    for _, row in df_results.iterrows():
        expected_total = row['n_outside_operational_scope'] + row['n_operationally_evaluable_but_not_assigned']
        actual_total = row['n_total_missing']
        if expected_total != actual_total:
            print(f"  ⚠ {row['criterion']}: Total mismatch! {expected_total} != {actual_total}")
    
    # Verify fractions sum to 1.0 for criteria with missing instances
    df_with_missing = df_results[df_results['n_total_missing'] > 0]
    for _, row in df_with_missing.iterrows():
        frac_sum = row['outside_scope_fraction_of_missing'] + row['operationally_evaluable_but_not_assigned_fraction']
        if not np.isclose(frac_sum, 1.0):
            print(f"  ⚠ {row['criterion']}: Fractions don't sum to 1.0! {frac_sum}")
    
    print("  ✓ All totals verified")
    print("  ✓ All fractions verified")
    
    # Save
    df_results.to_csv(DECOMP_FILE, index=False)
    print(f"\n✓ Saved: {DECOMP_FILE.name}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    n_with_missing = (df_results['n_total_missing'] > 0).sum()
    n_without_missing = (df_results['n_total_missing'] == 0).sum()
    
    print(f"\nCriteria with missing instances: {n_with_missing}")
    print(f"Criteria without missing instances: {n_without_missing}")
    
    # Top criteria by total missing
    print("\nTop 10 criteria by total missing:")
    df_top = df_results.nlargest(10, 'n_total_missing')
    print(f"\n{'Criterion':<10} {'Total':>8} {'Outside':>8} {'Evaluable':>10} {'Outside %':>10} {'Evaluable %':>12}")
    print("-" * 70)
    for _, row in df_top.iterrows():
        if row['n_total_missing'] > 0:
            print(f"{row['criterion']:<10} {row['n_total_missing']:>8} "
                  f"{row['n_outside_operational_scope']:>8} "
                  f"{row['n_operationally_evaluable_but_not_assigned']:>10} "
                  f"{row['outside_scope_fraction_of_missing']*100:>9.1f}% "
                  f"{row['operationally_evaluable_but_not_assigned_fraction']*100:>11.1f}%")
    
    # Criteria with 100% outside scope
    df_100_outside = df_results[
        (df_results['n_total_missing'] > 0) &
        (df_results['outside_scope_fraction_of_missing'] == 1.0)
    ]
    if len(df_100_outside) > 0:
        print(f"\nCriteria with 100% outside operational scope: {len(df_100_outside)}")
        for _, row in df_100_outside.iterrows():
            print(f"  {row['criterion']}: n={row['n_total_missing']}")
    
    # Criteria with 100% evaluable but not assigned
    df_100_evaluable = df_results[
        (df_results['n_total_missing'] > 0) &
        (df_results['operationally_evaluable_but_not_assigned_fraction'] == 1.0)
    ]
    if len(df_100_evaluable) > 0:
        print(f"\nCriteria with 100% operationally evaluable but not assigned: {len(df_100_evaluable)}")
        for _, row in df_100_evaluable.iterrows():
            print(f"  {row['criterion']}: n={row['n_total_missing']}")
    
    print("\n" + "="*80)
    print("✓ FOXL2 CRITERION DECOMPOSITION COMPLETE")
    print("="*80)
    print()
    
    return df_results


def main():
    """Main execution"""
    df_results = generate_criterion_decomposition()
    
    print(f"Output: {DECOMP_FILE}")
    print()


if __name__ == '__main__':
    main()
