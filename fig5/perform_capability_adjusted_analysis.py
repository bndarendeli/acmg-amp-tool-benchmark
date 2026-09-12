#!/usr/bin/env python3
"""
Capability-Adjusted Evidence Analysis for Figure 5

Performs capability-adjusted evidence concordance analysis using operational
criterion support from the capability matrix.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
import sys

# Add path for criteria parser
sys.path.insert(0, str(Path(__file__).parent.parent / 'fig1'))
from utils_criteria_parser import parse_criteria_string

# Configuration
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'outputs'

# Input files
CAPABILITY_MATRIX = OUTPUT_DIR / 'tool_criterion_capability_matrix.csv'
MASTER_TABLE = OUTPUT_DIR / 'fig5_master_variant_tool.csv'

# Output files
VARIANT_LEVEL_FILE = OUTPUT_DIR / 'fig5_capability_adjusted_variant_tool.csv'
SUMMARY_FILE = OUTPUT_DIR / 'fig5_capability_adjusted_summary.csv'
DECOMPOSITION_FILE = OUTPUT_DIR / 'fig5_missing_criteria_decomposition.csv'
REPORT_FILE = OUTPUT_DIR / 'fig5_capability_adjusted_report.txt'

# Canonical ACMG criteria
CANONICAL_CRITERIA = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]


def load_capability_matrix():
    """
    Load the capability matrix and create operational criterion sets per tool.
    
    Returns:
        dict: {tool: set of operationally supported criteria}
    """
    print("\nLoading capability matrix...")
    df_cap = pd.read_csv(CAPABILITY_MATRIX)
    print(f"✓ Loaded {len(df_cap)} tool×criterion pairs")
    
    # Build operational criterion sets
    operational_sets = {}
    for tool in df_cap['tool'].unique():
        df_tool = df_cap[df_cap['tool'] == tool]
        supported = set(df_tool[df_tool['operational_status'] == 'supported']['criterion'])
        operational_sets[tool] = supported
        print(f"  {tool}: {len(supported)}/28 operational criteria")
    
    return operational_sets


def calculate_capability_adjusted_metrics(G_i, T_it, S_t):
    """
    Calculate capability-adjusted metrics for a single variant-tool pair.
    
    Args:
        G_i: Reference criterion set (set or None)
        T_it: Tool criterion set (set or None)
        S_t: Operational criterion set for tool (set)
    
    Returns:
        dict: Capability-adjusted metrics
    """
    # Treat None as empty set (consistent with master table Jaccard calculation)
    if G_i is None:
        G_i = set()
    if T_it is None:
        T_it = set()
    
    # 1. Operationally eligible reference evidence
    G_supported = G_i & S_t
    
    # 2. Operationally eligible tool evidence
    T_supported = T_it & S_t
    
    # 3. Capability-adjusted Jaccard
    if len(G_supported) == 0 and len(T_supported) == 0:
        # Both empty - flag as no operationally evaluable criteria
        adjusted_jaccard = None
        flag = 'no_operationally_evaluable_criteria'
    elif len(G_supported) > 0 or len(T_supported) > 0:
        # At least one non-empty
        union = G_supported | T_supported
        intersection = G_supported & T_supported
        adjusted_jaccard = len(intersection) / len(union) if len(union) > 0 else 0.0
        flag = 'evaluable'
    else:
        adjusted_jaccard = None
        flag = 'evaluable'
    
    # 4. Capability-adjusted exact match
    if len(G_supported) == 0 and len(T_supported) == 0:
        adjusted_exact_match = None
    else:
        adjusted_exact_match = (G_supported == T_supported)
    
    # 5. Supported-reference recall
    if len(G_supported) > 0:
        intersection = G_i & T_it & S_t
        supported_reference_recall = len(intersection) / len(G_supported)
    else:
        supported_reference_recall = None
    
    # 6. Unsupported reference evidence fraction
    if len(G_i) > 0:
        unsupported = G_i - S_t
        unsupported_reference_fraction = len(unsupported) / len(G_i)
    else:
        unsupported_reference_fraction = None
    
    # 7. Decompose missing criteria
    missing = G_i - T_it
    
    # A. Outside operational scope
    outside_scope_missing = missing - S_t
    
    # B. Supported but missed
    supported_but_missed = missing & S_t
    
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
    """
    Perform capability-adjusted analysis at variant level.
    
    Args:
        df_master: Master variant-tool table
        operational_sets: dict of operational criterion sets per tool
    
    Returns:
        pd.DataFrame: Variant-level results
    """
    print("\nPerforming variant-level capability-adjusted analysis...")
    
    # Filter to ClinGen, classification-concordant, evidence-evaluable
    # IMPORTANT: Use SAME filter as original analysis (criterion_jaccard.notna())
    # to ensure population consistency
    df_clingen = df_master[df_master['dataset'] == 'clingen_28012026'].copy()
    print(f"  ClinGen predictions: {len(df_clingen):,}")
    
    df_concordant = df_clingen[df_clingen['classification_match'] == True].copy()
    print(f"  Classification-concordant: {len(df_concordant):,}")
    
    # Evidence-evaluable: criterion_jaccard is not NA
    # This matches the original evidence distribution analysis
    df_evaluable = df_concordant[df_concordant['criterion_jaccard'].notna()].copy()
    print(f"  Evidence-evaluable (criterion_jaccard notna): {len(df_evaluable):,}")
    
    # Perform capability-adjusted analysis
    results = []
    
    for idx, row in df_evaluable.iterrows():
        tool = row['tool']
        
        # Parse criteria (treat NA as empty set, consistent with master table)
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
            'flag': metrics['flag']
        }
        
        results.append(result)
    
    df_results = pd.DataFrame(results)
    
    print(f"✓ Analyzed {len(df_results):,} variant-tool pairs")
    
    # Report flags
    flag_counts = df_results['flag'].value_counts()
    print(f"\n  Flags:")
    for flag, count in flag_counts.items():
        print(f"    {flag}: {count:,}")
    
    return df_results


def calculate_tool_summary(df_variant):
    """
    Calculate tool-level summary statistics.
    
    Args:
        df_variant: Variant-level results
    
    Returns:
        pd.DataFrame: Tool summary
    """
    print("\nCalculating tool-level summary...")
    
    summaries = []
    
    for tool in sorted(df_variant['tool'].unique()):
        df_tool = df_variant[df_variant['tool'] == tool]
        
        # Total evaluable
        n_evaluable = len(df_tool)
        
        # Adjusted evaluable (exclude no_operationally_evaluable_criteria)
        df_adjusted_evaluable = df_tool[df_tool['flag'] == 'evaluable']
        n_adjusted_evaluable = len(df_adjusted_evaluable)
        
        # Raw metrics (all evaluable)
        raw_exact_match_rate = df_tool['raw_exact_match'].mean()
        raw_mean_jaccard = df_tool['raw_jaccard'].mean()
        raw_median_jaccard = df_tool['raw_jaccard'].median()
        
        # Adjusted metrics (only adjusted evaluable)
        if n_adjusted_evaluable > 0:
            adjusted_exact_match_rate = df_adjusted_evaluable['adjusted_exact_match'].mean()
            adjusted_mean_jaccard = df_adjusted_evaluable['adjusted_jaccard'].mean()
            adjusted_median_jaccard = df_adjusted_evaluable['adjusted_jaccard'].median()
            median_supported_reference_recall = df_adjusted_evaluable['supported_reference_recall'].median()
            mean_supported_reference_recall = df_adjusted_evaluable['supported_reference_recall'].mean()
            median_unsupported_reference_fraction = df_adjusted_evaluable['unsupported_reference_fraction'].median()
            mean_unsupported_reference_fraction = df_adjusted_evaluable['unsupported_reference_fraction'].mean()
        else:
            adjusted_exact_match_rate = None
            adjusted_mean_jaccard = None
            adjusted_median_jaccard = None
            median_supported_reference_recall = None
            mean_supported_reference_recall = None
            median_unsupported_reference_fraction = None
            mean_unsupported_reference_fraction = None
        
        # Deltas
        if adjusted_median_jaccard is not None and raw_median_jaccard is not None:
            delta_median_jaccard = adjusted_median_jaccard - raw_median_jaccard
        else:
            delta_median_jaccard = None
        
        if adjusted_exact_match_rate is not None and raw_exact_match_rate is not None:
            delta_exact_match_rate = adjusted_exact_match_rate - raw_exact_match_rate
        else:
            delta_exact_match_rate = None
        
        summaries.append({
            'tool': tool,
            'n_evaluable': n_evaluable,
            'n_adjusted_evaluable': n_adjusted_evaluable,
            'raw_exact_match_rate': raw_exact_match_rate,
            'adjusted_exact_match_rate': adjusted_exact_match_rate,
            'raw_mean_jaccard': raw_mean_jaccard,
            'raw_median_jaccard': raw_median_jaccard,
            'adjusted_mean_jaccard': adjusted_mean_jaccard,
            'adjusted_median_jaccard': adjusted_median_jaccard,
            'median_supported_reference_recall': median_supported_reference_recall,
            'mean_supported_reference_recall': mean_supported_reference_recall,
            'median_unsupported_reference_fraction': median_unsupported_reference_fraction,
            'mean_unsupported_reference_fraction': mean_unsupported_reference_fraction,
            'delta_median_jaccard': delta_median_jaccard,
            'delta_exact_match_rate': delta_exact_match_rate
        })
    
    df_summary = pd.DataFrame(summaries)
    
    print(f"✓ Summarized {len(df_summary)} tools")
    
    return df_summary


def decompose_missing_criteria(df_variant):
    """
    Decompose missing criteria into outside-scope vs. supported-but-missed.
    
    Args:
        df_variant: Variant-level results
    
    Returns:
        tuple: (pooled_df, tool_specific_df)
    """
    print("\nDecomposing missing criteria...")
    
    # Pooled analysis
    pooled_results = []
    
    for criterion in CANONICAL_CRITERIA:
        # Count across all variant-tool pairs
        n_outside_scope = 0
        n_supported_but_missed = 0
        
        for _, row in df_variant.iterrows():
            outside = row['outside_scope_missing_criteria']
            supported_missed = row['supported_but_missed_criteria']
            
            if pd.notna(outside) and outside:
                if criterion in outside.split(','):
                    n_outside_scope += 1
            
            if pd.notna(supported_missed) and supported_missed:
                if criterion in supported_missed.split(','):
                    n_supported_but_missed += 1
        
        n_total_missing = n_outside_scope + n_supported_but_missed
        
        if n_total_missing > 0:
            outside_scope_fraction = n_outside_scope / n_total_missing
            supported_but_missed_fraction = n_supported_but_missed / n_total_missing
        else:
            outside_scope_fraction = 0.0
            supported_but_missed_fraction = 0.0
        
        pooled_results.append({
            'criterion': criterion,
            'n_total_missing': n_total_missing,
            'n_outside_operational_scope': n_outside_scope,
            'n_supported_but_missed': n_supported_but_missed,
            'outside_scope_fraction_of_missing': outside_scope_fraction,
            'supported_but_missed_fraction_of_missing': supported_but_missed_fraction
        })
    
    df_pooled = pd.DataFrame(pooled_results)
    df_pooled = df_pooled.sort_values('n_total_missing', ascending=False)
    
    # Tool-specific analysis
    tool_specific_results = []
    
    for tool in sorted(df_variant['tool'].unique()):
        df_tool = df_variant[df_variant['tool'] == tool]
        
        for criterion in CANONICAL_CRITERIA:
            n_outside_scope = 0
            n_supported_but_missed = 0
            
            for _, row in df_tool.iterrows():
                outside = row['outside_scope_missing_criteria']
                supported_missed = row['supported_but_missed_criteria']
                
                if pd.notna(outside) and outside:
                    if criterion in outside.split(','):
                        n_outside_scope += 1
                
                if pd.notna(supported_missed) and supported_missed:
                    if criterion in supported_missed.split(','):
                        n_supported_but_missed += 1
            
            n_total_missing = n_outside_scope + n_supported_but_missed
            
            if n_total_missing > 0:
                outside_scope_fraction = n_outside_scope / n_total_missing
                supported_but_missed_fraction = n_supported_but_missed / n_total_missing
            else:
                outside_scope_fraction = 0.0
                supported_but_missed_fraction = 0.0
            
            tool_specific_results.append({
                'tool': tool,
                'criterion': criterion,
                'n_total_missing': n_total_missing,
                'n_outside_operational_scope': n_outside_scope,
                'n_supported_but_missed': n_supported_but_missed,
                'outside_scope_fraction_of_missing': outside_scope_fraction,
                'supported_but_missed_fraction_of_missing': supported_but_missed_fraction
            })
    
    df_tool_specific = pd.DataFrame(tool_specific_results)
    
    print(f"✓ Decomposed {len(df_pooled)} criteria (pooled)")
    print(f"✓ Decomposed {len(df_tool_specific)} tool×criterion pairs")
    
    return df_pooled, df_tool_specific


def generate_report(df_summary, df_pooled, df_tool_specific):
    """
    Generate capability-adjusted analysis report.
    
    Args:
        df_summary: Tool summary
        df_pooled: Pooled criterion decomposition
        df_tool_specific: Tool-specific criterion decomposition
    
    Returns:
        str: Report text
    """
    lines = []
    lines.append("="*80)
    lines.append("CAPABILITY-ADJUSTED EVIDENCE ANALYSIS REPORT")
    lines.append("="*80)
    lines.append("")
    
    # Section 1: Overview
    lines.append("SECTION 1: ANALYSIS OVERVIEW")
    lines.append("-"*80)
    lines.append("")
    lines.append("Population:")
    lines.append("  Dataset: ClinGen (clingen_28012026)")
    lines.append("  Filter: Classification-concordant, evidence-evaluable predictions")
    lines.append("")
    lines.append("Capability Adjustment:")
    lines.append("  Used operational_status = 'supported' from capability matrix")
    lines.append("  Restricted evidence comparison to operationally supported criteria")
    lines.append("")
    
    # Section 2: Raw vs. Adjusted Comparison
    lines.append("")
    lines.append("SECTION 2: RAW VS. CAPABILITY-ADJUSTED EVIDENCE FIDELITY")
    lines.append("-"*80)
    lines.append("")
    
    # Sort by delta_median_jaccard
    df_sorted = df_summary.sort_values('delta_median_jaccard', ascending=False)
    
    lines.append(f"{'Tool':<20s} {'Raw J':<8s} {'Adj J':<8s} {'ΔJ':<8s} {'Raw Exact':<10s} {'Adj Exact':<10s} {'ΔExact':<8s}")
    lines.append("-"*80)
    
    for _, row in df_sorted.iterrows():
        tool = row['tool']
        raw_j = row['raw_median_jaccard']
        adj_j = row['adjusted_median_jaccard']
        delta_j = row['delta_median_jaccard']
        raw_exact = row['raw_exact_match_rate']
        adj_exact = row['adjusted_exact_match_rate']
        delta_exact = row['delta_exact_match_rate']
        
        raw_j_str = f"{raw_j:.3f}" if pd.notna(raw_j) else "N/A"
        adj_j_str = f"{adj_j:.3f}" if pd.notna(adj_j) else "N/A"
        delta_j_str = f"{delta_j:+.3f}" if pd.notna(delta_j) else "N/A"
        raw_exact_str = f"{raw_exact*100:.1f}%" if pd.notna(raw_exact) else "N/A"
        adj_exact_str = f"{adj_exact*100:.1f}%" if pd.notna(adj_exact) else "N/A"
        delta_exact_str = f"{delta_exact*100:+.1f}%" if pd.notna(delta_exact) else "N/A"
        
        lines.append(f"{tool:<20s} {raw_j_str:<8s} {adj_j_str:<8s} {delta_j_str:<8s} {raw_exact_str:<10s} {adj_exact_str:<10s} {delta_exact_str:<8s}")
    
    lines.append("")
    lines.append("Legend:")
    lines.append("  Raw J: Raw median Jaccard similarity")
    lines.append("  Adj J: Capability-adjusted median Jaccard similarity")
    lines.append("  ΔJ: Change in median Jaccard (adjusted - raw)")
    lines.append("  Raw Exact: Raw exact evidence match rate")
    lines.append("  Adj Exact: Capability-adjusted exact evidence match rate")
    lines.append("  ΔExact: Change in exact match rate (adjusted - raw)")
    lines.append("")
    
    # Section 3: Tool-Specific Analysis
    lines.append("")
    lines.append("SECTION 3: TOOL-SPECIFIC CAPABILITY-ADJUSTED ANALYSIS")
    lines.append("-"*80)
    lines.append("")
    
    for _, row in df_summary.sort_values('tool').iterrows():
        tool = row['tool']
        lines.append(f"{tool}:")
        lines.append(f"  Evidence-evaluable predictions: {row['n_evaluable']:,}")
        lines.append(f"  Adjusted-evaluable predictions: {row['n_adjusted_evaluable']:,}")
        lines.append("")
        adj_j_val = f"{row['adjusted_median_jaccard']:.3f}" if pd.notna(row['adjusted_median_jaccard']) else 'N/A'
        delta_j_val = f"{row['delta_median_jaccard']:+.3f}" if pd.notna(row['delta_median_jaccard']) else 'N/A'
        recall_val = f"{row['median_supported_reference_recall']:.3f}" if pd.notna(row['median_supported_reference_recall']) else 'N/A'
        unsup_val = f"{row['median_unsupported_reference_fraction']:.3f}" if pd.notna(row['median_unsupported_reference_fraction']) else 'N/A'
        
        lines.append(f"  Raw median Jaccard: {row['raw_median_jaccard']:.3f}")
        lines.append(f"  Adjusted median Jaccard: {adj_j_val}")
        lines.append(f"  Change: {delta_j_val}")
        lines.append("")
        lines.append(f"  Median supported-reference recall: {recall_val}")
        lines.append(f"  Median unsupported-reference fraction: {unsup_val}")
        lines.append("")
        
        # Characterize low fidelity
        if pd.notna(row['median_supported_reference_recall']) and pd.notna(row['median_unsupported_reference_fraction']):
            recall = row['median_supported_reference_recall']
            unsupported_frac = row['median_unsupported_reference_fraction']
            
            if recall < 0.5 and unsupported_frac > 0.3:
                characterization = "Both limited operational scope AND within-scope disagreement"
            elif unsupported_frac > 0.3:
                characterization = "Predominantly limited operational criterion scope"
            elif recall < 0.5:
                characterization = "Predominantly within-scope assignment disagreement"
            else:
                characterization = "Good operational coverage and within-scope agreement"
            
            lines.append(f"  Characterization: {characterization}")
        
        lines.append("")
    
    # Section 4: Missing Criteria Decomposition
    lines.append("")
    lines.append("SECTION 4: MISSING CRITERIA DECOMPOSITION (POOLED)")
    lines.append("-"*80)
    lines.append("")
    lines.append("Top 10 most frequently missing criteria:")
    lines.append("")
    lines.append(f"{'Criterion':<12s} {'Total':<8s} {'Outside Scope':<15s} {'Supported/Missed':<15s} {'% Outside':<12s}")
    lines.append("-"*80)
    
    for _, row in df_pooled.head(10).iterrows():
        criterion = row['criterion']
        total = row['n_total_missing']
        outside = row['n_outside_operational_scope']
        supported = row['n_supported_but_missed']
        pct_outside = row['outside_scope_fraction_of_missing'] * 100
        
        lines.append(f"{criterion:<12s} {total:<8d} {outside:<15d} {supported:<15d} {pct_outside:>10.1f}%")
    
    lines.append("")
    lines.append("Interpretation:")
    lines.append("  Outside Scope: Criterion not operationally supported by tool")
    lines.append("  Supported/Missed: Criterion operationally supported but not assigned")
    lines.append("")
    
    # Highlight major criteria
    lines.append("")
    lines.append("Previously Identified Major Missing Criteria:")
    lines.append("")
    
    major_criteria = ['PP4', 'PS4', 'PM3', 'BP4', 'PP1']
    for criterion in major_criteria:
        row = df_pooled[df_pooled['criterion'] == criterion].iloc[0]
        total = row['n_total_missing']
        outside = row['n_outside_operational_scope']
        supported = row['n_supported_but_missed']
        pct_outside = row['outside_scope_fraction_of_missing'] * 100
        
        lines.append(f"  {criterion}:")
        lines.append(f"    Total missing: {total:,}")
        lines.append(f"    Outside operational scope: {outside:,} ({pct_outside:.1f}%)")
        lines.append(f"    Supported but not assigned: {supported:,} ({100-pct_outside:.1f}%)")
        lines.append("")
    
    # Section 5: Key Findings
    lines.append("")
    lines.append("SECTION 5: KEY FINDINGS")
    lines.append("-"*80)
    lines.append("")
    
    # Calculate overall statistics
    mean_delta_j = df_summary['delta_median_jaccard'].mean()
    mean_delta_exact = df_summary['delta_exact_match_rate'].mean()
    
    lines.append("1. Overall Impact of Capability Adjustment:")
    lines.append(f"   Mean change in median Jaccard: {mean_delta_j:+.3f}")
    lines.append(f"   Mean change in exact match rate: {mean_delta_exact*100:+.1f}%")
    lines.append("")
    
    # Tools with largest improvement
    top_improved = df_summary.nlargest(3, 'delta_median_jaccard')
    lines.append("2. Tools with Largest Improvement (capability adjustment):")
    for _, row in top_improved.iterrows():
        lines.append(f"   {row['tool']}: ΔJ = {row['delta_median_jaccard']:+.3f}")
    lines.append("")
    
    # Tools with smallest change
    smallest_change = df_summary.nsmallest(3, 'delta_median_jaccard', keep='all')
    lines.append("3. Tools with Smallest Change (already well-matched to operational scope):")
    for _, row in smallest_change.iterrows():
        lines.append(f"   {row['tool']}: ΔJ = {row['delta_median_jaccard']:+.3f}")
    lines.append("")
    
    lines.append("="*80)
    lines.append("END OF REPORT")
    lines.append("="*80)
    
    return '\n'.join(lines)


def main():
    """Main execution"""
    print("="*80)
    print("CAPABILITY-ADJUSTED EVIDENCE ANALYSIS FOR FIGURE 4")
    print("="*80)
    
    # Load capability matrix
    operational_sets = load_capability_matrix()
    
    # Load master table
    print(f"\nLoading master table: {MASTER_TABLE}")
    df_master = pd.read_csv(MASTER_TABLE)
    print(f"✓ Loaded {len(df_master):,} rows")
    
    # Perform variant-level analysis
    df_variant = perform_variant_level_analysis(df_master, operational_sets)
    
    # Calculate tool summary
    df_summary = calculate_tool_summary(df_variant)
    
    # Decompose missing criteria
    df_pooled, df_tool_specific = decompose_missing_criteria(df_variant)
    
    # Combine decomposition outputs
    df_decomposition = pd.concat([
        df_pooled.assign(analysis_type='pooled'),
        df_tool_specific.assign(analysis_type='tool_specific')
    ], ignore_index=True)
    
    # Generate report
    print("\nGenerating report...")
    report = generate_report(df_summary, df_pooled, df_tool_specific)
    
    # Save outputs
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\nSaving outputs...")
    
    print(f"  Variant-level: {VARIANT_LEVEL_FILE}")
    df_variant.to_csv(VARIANT_LEVEL_FILE, index=False)
    print(f"    ✓ Saved {len(df_variant):,} rows")
    
    print(f"  Tool summary: {SUMMARY_FILE}")
    df_summary.to_csv(SUMMARY_FILE, index=False)
    print(f"    ✓ Saved {len(df_summary)} tools")
    
    print(f"  Criterion decomposition: {DECOMPOSITION_FILE}")
    df_decomposition.to_csv(DECOMPOSITION_FILE, index=False)
    print(f"    ✓ Saved {len(df_decomposition)} rows")
    
    print(f"  Report: {REPORT_FILE}")
    with open(REPORT_FILE, 'w') as f:
        f.write(report)
    print(f"    ✓ Saved report")
    
    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    print(f"\nCapability-Adjusted Analysis Complete:")
    print(f"  Variant-tool pairs analyzed: {len(df_variant):,}")
    print(f"  Tools summarized: {len(df_summary)}")
    print(f"  Criteria decomposed: {len(df_pooled)}")
    
    print(f"\nOverall Impact:")
    mean_delta_j = df_summary['delta_median_jaccard'].mean()
    mean_delta_exact = df_summary['delta_exact_match_rate'].mean()
    print(f"  Mean ΔJaccard: {mean_delta_j:+.3f}")
    print(f"  Mean ΔExact: {mean_delta_exact*100:+.1f}%")
    
    print("\n" + "="*80)
    print("✓ CAPABILITY-ADJUSTED ANALYSIS COMPLETE")
    print("="*80)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
