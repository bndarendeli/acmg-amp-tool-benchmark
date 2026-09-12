#!/usr/bin/env python3
"""
FOXL2 Capability-Adjusted Evidence Analysis

Performs capability-adjusted evidence concordance analysis for FOXL2 dataset
using the same methodology as the corrected ClinGen analysis.
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
MASTER_TABLE = BASE_DIR / 'outputs' / 'fig5_master_variant_tool.csv'

# Output files
VARIANT_LEVEL_FILE = OUTPUT_DIR / 'foxl2_capability_adjusted_variant_tool.csv'
SUMMARY_FILE = OUTPUT_DIR / 'foxl2_capability_adjusted_summary.csv'


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


def calculate_capability_adjusted_metrics(G_i, T_it, S_t):
    """
    Calculate capability-adjusted metrics for a single variant-tool pair.
    
    Args:
        G_i: Reference criterion set
        T_it: Tool criterion set
        S_t: Tool's operational criterion set
    
    Returns:
        dict: Adjusted metrics
    """
    # Project onto operational scope
    G_supported = G_i.intersection(S_t) if G_i is not None else set()
    T_supported = T_it.intersection(S_t) if T_it is not None else set()
    
    # Calculate adjusted Jaccard
    if len(G_supported) == 0 and len(T_supported) == 0:
        # Both projected sets empty
        flag = 'no_operationally_evaluable_criteria'
        adjusted_jaccard = None
        adjusted_exact_match = None
    else:
        flag = 'evaluable'
        # Standard Jaccard
        union = G_supported.union(T_supported)
        if len(union) == 0:
            adjusted_jaccard = 1.0
        else:
            intersection = G_supported.intersection(T_supported)
            adjusted_jaccard = len(intersection) / len(union)
        
        adjusted_exact_match = (G_supported == T_supported)
    
    # Decompose missing criteria
    missing = G_i.difference(T_it) if G_i is not None and T_it is not None else set()
    outside_scope_missing = missing.difference(S_t)
    supported_but_missed = missing.intersection(S_t)
    
    # Calculate supported-reference recall: |G ∩ T ∩ S| / |G ∩ S|
    # This measures how well the tool reproduces the supported reference criteria
    if len(G_supported) > 0:
        # Intersection of reference and tool within operational scope
        supported_intersection = G_supported.intersection(T_supported)
        supported_reference_recall = len(supported_intersection) / len(G_supported)
    else:
        supported_reference_recall = None
    
    # Calculate unsupported reference fraction for context
    if G_i is not None and len(G_i) > 0:
        unsupported_reference_fraction = len(G_i.difference(S_t)) / len(G_i)
    else:
        unsupported_reference_fraction = None
    
    return {
        'supported_reference_criteria': sorted(G_supported),
        'supported_tool_criteria': sorted(T_supported),
        'adjusted_jaccard': adjusted_jaccard,
        'adjusted_exact_match': adjusted_exact_match,
        'supported_reference_recall': supported_reference_recall,
        'unsupported_reference_fraction': unsupported_reference_fraction,
        'outside_scope_missing': sorted(outside_scope_missing),
        'supported_but_missed': sorted(supported_but_missed),
        'n_outside_scope_missing': len(outside_scope_missing),
        'n_supported_but_missed': len(supported_but_missed),
        'flag': flag
    }


def perform_variant_level_analysis(df_master, operational_sets):
    """Perform capability-adjusted analysis at variant level for FOXL2"""
    print("\nPerforming FOXL2 variant-level capability-adjusted analysis...")
    
    # Filter to FOXL2, classification-concordant, evidence-evaluable
    df_foxl2 = df_master[df_master['dataset'] == 'foxl2'].copy()
    print(f"  FOXL2 predictions: {len(df_foxl2):,}")
    
    df_concordant = df_foxl2[df_foxl2['classification_match'] == True].copy()
    print(f"  Classification-concordant: {len(df_concordant):,}")
    
    # Evidence-evaluable: criterion_jaccard is not NA
    df_evaluable = df_concordant[df_concordant['criterion_jaccard'].notna()].copy()
    print(f"  Evidence-evaluable (criterion_jaccard notna): {len(df_evaluable):,}")
    
    # Perform capability-adjusted analysis
    results = []
    
    for idx, row in df_evaluable.iterrows():
        tool = row['tool']
        
        # Parse criteria (treat NA as empty set)
        ref_crit = row['reference_criteria']
        tool_crit = row['tool_criteria']
        
        if pd.isna(ref_crit):
            G_i = set()
        else:
            G_i = parse_criteria_string(ref_crit)
        
        if pd.isna(tool_crit):
            T_it = set()
        else:
            T_it = parse_criteria_string(tool_crit)
        
        S_t = operational_sets.get(tool, set())
        
        # Calculate adjusted metrics
        metrics = calculate_capability_adjusted_metrics(G_i, T_it, S_t)
        
        # Build result row
        result = {
            'dataset': row['dataset'],
            'variant_id': row['variant_id'],
            'tool': tool,
            'reference_class': row['reference_class'],
            'tool_class': row['tool_class'],
            'raw_reference_criteria': row['reference_criteria'],
            'raw_tool_criteria': row['tool_criteria'],
            'operational_supported_criteria': ','.join(sorted(S_t)),
            'supported_reference_criteria': ','.join(metrics['supported_reference_criteria']) if metrics['supported_reference_criteria'] else '',
            'supported_tool_criteria': ','.join(metrics['supported_tool_criteria']) if metrics['supported_tool_criteria'] else '',
            'raw_jaccard': row['criterion_jaccard'],
            'adjusted_jaccard': metrics['adjusted_jaccard'],
            'raw_exact_match': row['exact_evidence_match'],
            'adjusted_exact_match': metrics['adjusted_exact_match'],
            'supported_reference_recall': metrics['supported_reference_recall'],
            'unsupported_reference_fraction': metrics['unsupported_reference_fraction'],
            'outside_scope_missing_criteria': ','.join(metrics['outside_scope_missing']) if metrics['outside_scope_missing'] else '',
            'supported_but_missed_criteria': ','.join(metrics['supported_but_missed']) if metrics['supported_but_missed'] else '',
            'n_outside_scope_missing': metrics['n_outside_scope_missing'],
            'n_supported_but_missed': metrics['n_supported_but_missed'],
            'no_operationally_evaluable_criteria': (metrics['flag'] == 'no_operationally_evaluable_criteria')
        }
        
        results.append(result)
    
    df_results = pd.DataFrame(results)
    
    print(f"✓ Analyzed {len(df_results):,} variant-tool pairs")
    
    # Report flags
    n_no_op = (df_results['no_operationally_evaluable_criteria'] == True).sum()
    n_evaluable = (df_results['no_operationally_evaluable_criteria'] == False).sum()
    print(f"\n  no_operationally_evaluable_criteria: {n_no_op:,}")
    print(f"  evaluable: {n_evaluable:,}")
    
    return df_results


def calculate_tool_summary(df_variant):
    """Calculate tool-level summary statistics"""
    print("\nCalculating tool-level summary...")
    
    tools = df_variant['tool'].unique()
    summaries = []
    
    for tool in tools:
        df_tool = df_variant[df_variant['tool'] == tool]
        
        # Filter to evaluable (not flagged as no_operationally_evaluable_criteria)
        df_evaluable = df_tool[df_tool['no_operationally_evaluable_criteria'] == False]
        
        n_total = len(df_tool)
        n_evaluable = len(df_evaluable)
        n_no_op = len(df_tool[df_tool['no_operationally_evaluable_criteria'] == True])
        
        if n_evaluable > 0:
            # Raw metrics
            raw_exact_match_rate = (df_evaluable['raw_exact_match'] == True).sum() / n_evaluable
            raw_mean_jaccard = df_evaluable['raw_jaccard'].mean()
            raw_median_jaccard = df_evaluable['raw_jaccard'].median()
            
            # Adjusted metrics
            adjusted_exact_match_rate = (df_evaluable['adjusted_exact_match'] == True).sum() / n_evaluable
            adjusted_mean_jaccard = df_evaluable['adjusted_jaccard'].mean()
            adjusted_median_jaccard = df_evaluable['adjusted_jaccard'].median()
            
            # Deltas
            delta_exact_match = adjusted_exact_match_rate - raw_exact_match_rate
            delta_mean_jaccard = adjusted_mean_jaccard - raw_mean_jaccard
            delta_median_jaccard = adjusted_median_jaccard - raw_median_jaccard
            
            # Recall metrics
            median_supported_recall = df_evaluable['supported_reference_recall'].median()
            mean_supported_recall = df_evaluable['supported_reference_recall'].mean()
            median_unsupported_frac = df_evaluable['unsupported_reference_fraction'].median()
            mean_unsupported_frac = df_evaluable['unsupported_reference_fraction'].mean()
        else:
            raw_exact_match_rate = np.nan
            raw_mean_jaccard = np.nan
            raw_median_jaccard = np.nan
            adjusted_exact_match_rate = np.nan
            adjusted_mean_jaccard = np.nan
            adjusted_median_jaccard = np.nan
            delta_exact_match = np.nan
            delta_mean_jaccard = np.nan
            delta_median_jaccard = np.nan
            median_supported_recall = np.nan
            mean_supported_recall = np.nan
            median_unsupported_frac = np.nan
            mean_unsupported_frac = np.nan
        
        summaries.append({
            'tool': tool,
            'n_evaluable': n_total,
            'n_adjusted_evaluable': n_evaluable,
            'n_no_operationally_evaluable_criteria': n_no_op,
            'raw_exact_match_rate': raw_exact_match_rate,
            'adjusted_exact_match_rate': adjusted_exact_match_rate,
            'raw_mean_jaccard': raw_mean_jaccard,
            'raw_median_jaccard': raw_median_jaccard,
            'adjusted_mean_jaccard': adjusted_mean_jaccard,
            'adjusted_median_jaccard': adjusted_median_jaccard,
            'median_supported_reference_recall': median_supported_recall,
            'mean_supported_reference_recall': mean_supported_recall,
            'median_unsupported_reference_fraction': median_unsupported_frac,
            'mean_unsupported_reference_fraction': mean_unsupported_frac,
            'delta_median_jaccard': delta_median_jaccard,
            'delta_mean_jaccard': delta_mean_jaccard,
            'delta_exact_match_rate': delta_exact_match
        })
    
    df_summary = pd.DataFrame(summaries)
    
    print(f"✓ Summarized {len(df_summary)} tools")
    
    return df_summary


def main():
    """Main execution"""
    print("="*80)
    print("FOXL2 CAPABILITY-ADJUSTED ANALYSIS")
    print("="*80)
    print()
    print("Using same methodology as corrected ClinGen analysis")
    print()
    
    # Load data
    df_master = pd.read_csv(MASTER_TABLE)
    operational_sets = load_capability_matrix()
    
    # Perform variant-level analysis
    df_variant = perform_variant_level_analysis(df_master, operational_sets)
    
    # Save variant-level results
    df_variant.to_csv(VARIANT_LEVEL_FILE, index=False)
    print(f"\n✓ Saved: {VARIANT_LEVEL_FILE.name}")
    
    # Calculate tool summary
    df_summary = calculate_tool_summary(df_variant)
    
    # Save summary
    df_summary.to_csv(SUMMARY_FILE, index=False)
    print(f"✓ Saved: {SUMMARY_FILE.name}")
    
    print()
    print("="*80)
    print("✓ FOXL2 CAPABILITY-ADJUSTED ANALYSIS COMPLETE")
    print("="*80)
    print()


if __name__ == '__main__':
    main()
