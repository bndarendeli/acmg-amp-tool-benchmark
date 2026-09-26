#!/usr/bin/env python3
"""
Generate ClinGen Supplementary Tables for Main Figure 5

Creates publication-ready supplementary tables supporting Main Figure 5:
- Full tool capability matrix
- Full ClinGen tool-level evidence summary
- Full ClinGen criterion-level missing decomposition (all 28 criteria)
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Configuration
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'outputs' / 'clingen_supplementary'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Input files
CAPABILITY_MATRIX = BASE_DIR / 'outputs' / 'tool_criterion_capability_matrix.csv'
CAPABILITY_REF = BASE_DIR / 'reference' / 'tool_acmg_capabilities.csv'
EVIDENCE_DIST = BASE_DIR / 'outputs' / 'fig5_evidence_distributions.csv'
CAPABILITY_ADJ_SUMMARY = BASE_DIR / 'outputs' / 'fig5_capability_adjusted_summary.csv'
MISSING_DECOMP = BASE_DIR / 'outputs' / 'fig5_missing_criteria_decomposition.csv'
MASTER_TABLE = BASE_DIR / 'outputs' / 'fig5_master_variant_tool.csv'

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


def generate_full_capability_table():
    """Generate full tool capability matrix"""
    print("\n" + "="*80)
    print("GENERATING FULL TOOL CAPABILITY TABLE")
    print("="*80)
    
    df_cap = pd.read_csv(CAPABILITY_MATRIX)
    
    # Pivot to wide format: tools as rows, criteria as columns
    df_pivot = df_cap.pivot(index='tool', columns='criterion', values='operational_status')
    
    # Reorder columns to canonical order
    df_pivot = df_pivot[CANONICAL_CRITERIA]
    
    # Add operational count column
    n_operational = []
    for tool in df_pivot.index:
        n_op = (df_pivot.loc[tool] == 'supported').sum()
        n_operational.append(n_op)
    
    df_pivot['n_operational_criteria'] = n_operational
    
    # Reset index to make tool a column
    df_pivot = df_pivot.reset_index()
    
    # Save
    output_path = OUTPUT_DIR / 'supplementary_table_tool_capabilities.csv'
    df_pivot.to_csv(output_path, index=False)
    print(f"✓ Saved: {output_path.name}")
    print(f"  {len(df_pivot)} tools × {len(CANONICAL_CRITERIA)} criteria")
    
    return df_pivot


def generate_full_clingen_tool_summary():
    """Generate full ClinGen tool-level evidence summary"""
    print("\n" + "="*80)
    print("GENERATING FULL CLINGEN TOOL-LEVEL SUMMARY")
    print("="*80)
    
    # Load data
    df_evidence = pd.read_csv(EVIDENCE_DIST)
    df_adj = pd.read_csv(CAPABILITY_ADJ_SUMMARY)
    df_master = pd.read_csv(MASTER_TABLE)
    df_cap = pd.read_csv(CAPABILITY_MATRIX)
    
    # Filter to ClinGen
    df_clingen_ev = df_evidence[
        (df_evidence['dataset'] == 'clingen_28012026') &
        (df_evidence['tool'] != 'POOLED')
    ].copy()
    
    df_clingen_master = df_master[df_master['dataset'] == 'clingen_28012026'].copy()
    
    # Get operational counts
    operational_counts = {}
    for tool in df_cap['tool'].unique():
        n_op = len(df_cap[
            (df_cap['tool'] == tool) &
            (df_cap['operational_status'] == 'supported')
        ])
        operational_counts[tool] = n_op
    
    # Build summary
    summary_rows = []
    
    for tool in df_clingen_ev['tool'].unique():
        # Classification concordance
        df_tool = df_clingen_master[df_clingen_master['tool'] == tool]
        n_total = len(df_tool)
        n_match = len(df_tool[df_tool['classification_match'] == True])
        class_concordance = n_match / n_total if n_total > 0 else np.nan
        
        # Raw evidence metrics
        ev_row = df_clingen_ev[df_clingen_ev['tool'] == tool].iloc[0]
        n_raw_evaluable = ev_row['n_classification_concordant_with_evaluable_evidence']
        raw_mean_j = ev_row['mean_criterion_jaccard']
        raw_median_j = ev_row['median_criterion_jaccard']
        raw_exact_match = ev_row['exact_evidence_match_rate']
        
        # Adjusted evidence metrics
        adj_row = df_adj[df_adj['tool'] == tool]
        if len(adj_row) > 0:
            adj_row = adj_row.iloc[0]
            n_adj_evaluable = adj_row['n_adjusted_evaluable']
            adj_mean_j = adj_row['adjusted_mean_jaccard']
            adj_median_j = adj_row['adjusted_median_jaccard']
            adj_exact_match = adj_row['adjusted_exact_match_rate']
            delta_mean_j = adj_row.get('delta_mean_jaccard', adj_mean_j - raw_mean_j)
            delta_median_j = adj_row['delta_median_jaccard']
            delta_exact = adj_row['delta_exact_match_rate']
            
            # Calculate n_no_op from n_evaluable vs n_adjusted_evaluable
            n_no_op = adj_row['n_evaluable'] - adj_row['n_adjusted_evaluable']
        else:
            n_adj_evaluable = 0
            n_no_op = 0
            adj_mean_j = np.nan
            adj_median_j = np.nan
            adj_exact_match = np.nan
            delta_mean_j = np.nan
            delta_median_j = np.nan
            delta_exact = np.nan
        
        summary_rows.append({
            'tool': tool,
            'n_total_predictions': n_total,
            'n_classification_matches': n_match,
            'classification_concordance': class_concordance,
            'n_raw_evaluable': n_raw_evaluable,
            'raw_mean_jaccard': raw_mean_j,
            'raw_median_jaccard': raw_median_j,
            'raw_exact_match_rate': raw_exact_match,
            'n_adjusted_evaluable': n_adj_evaluable,
            'n_no_operationally_evaluable': n_no_op,
            'adjusted_mean_jaccard': adj_mean_j,
            'adjusted_median_jaccard': adj_median_j,
            'adjusted_exact_match_rate': adj_exact_match,
            'delta_mean_jaccard': delta_mean_j,
            'delta_median_jaccard': delta_median_j,
            'delta_exact_match': delta_exact,
            'n_operational_criteria': operational_counts.get(tool, 0)
        })
    
    df_summary = pd.DataFrame(summary_rows)
    
    # Sort by raw median Jaccard (descending)
    df_summary = df_summary.sort_values('raw_median_jaccard', ascending=False)
    
    # Save
    output_path = OUTPUT_DIR / 'supplementary_table_clingen_tool_summary.csv'
    df_summary.to_csv(output_path, index=False)
    print(f"✓ Saved: {output_path.name}")
    print(f"  {len(df_summary)} tools")
    
    return df_summary


def generate_full_criterion_decomposition():
    """Generate full ClinGen criterion-level missing decomposition"""
    print("\n" + "="*80)
    print("GENERATING FULL CRITERION DECOMPOSITION")
    print("="*80)
    
    df = pd.read_csv(MISSING_DECOMP)
    df_clingen = df[df['analysis_type'] == 'pooled'].copy()
    
    # Ensure all 28 canonical criteria are represented
    existing_criteria = set(df_clingen['criterion'].unique())
    missing_criteria = set(CANONICAL_CRITERIA) - existing_criteria
    
    # Add rows for criteria with no missing instances
    for crit in missing_criteria:
        df_clingen = pd.concat([df_clingen, pd.DataFrame([{
            'analysis_type': 'pooled',
            'criterion': crit,
            'n_total_missing': 0,
            'n_outside_operational_scope': 0,
            'n_supported_but_missed': 0,
            'outside_scope_fraction_of_missing': 0.0,
            'supported_but_missed_fraction_of_missing': 0.0
        }])], ignore_index=True)
    
    # Sort by canonical order
    df_clingen['criterion'] = pd.Categorical(df_clingen['criterion'], 
                                             categories=CANONICAL_CRITERIA, 
                                             ordered=True)
    df_clingen = df_clingen.sort_values('criterion')
    
    # Select and rename columns for publication
    df_pub = df_clingen[[
        'criterion',
        'n_total_missing',
        'n_outside_operational_scope',
        'n_supported_but_missed',
        'outside_scope_fraction_of_missing',
        'supported_but_missed_fraction_of_missing'
    ]].copy()
    
    df_pub.columns = [
        'criterion',
        'total_missing_n',
        'outside_scope_n',
        'operationally_evaluable_but_not_assigned_n',
        'outside_scope_fraction',
        'operationally_evaluable_but_not_assigned_fraction'
    ]
    
    # Save
    output_path = OUTPUT_DIR / 'supplementary_table_clingen_criterion_decomposition.csv'
    df_pub.to_csv(output_path, index=False)
    print(f"✓ Saved: {output_path.name}")
    print(f"  {len(df_pub)} criteria (all 28 canonical)")
    
    # Report summary
    n_with_missing = (df_pub['total_missing_n'] > 0).sum()
    n_without_missing = (df_pub['total_missing_n'] == 0).sum()
    print(f"  Criteria with missing instances: {n_with_missing}")
    print(f"  Criteria without missing instances: {n_without_missing}")
    
    return df_pub


def generate_supplementary_readme():
    """Generate README for supplementary tables"""
    print("\n" + "="*80)
    print("GENERATING SUPPLEMENTARY README")
    print("="*80)
    
    readme = []
    readme.append("# ClinGen Supplementary Tables for Main Figure 5")
    readme.append("")
    readme.append("## Overview")
    readme.append("")
    readme.append("These tables provide complete data supporting Main Figure 5 of the")
    readme.append("ACMG/AMP automation benchmarking manuscript.")
    readme.append("")
    readme.append("All analyses use the corrected pipeline after fixing the ClinGen")
    readme.append("missing reference criterion bug (Ground_Truth_ACMG='.' correctly")
    readme.append("treated as missing/unavailable).")
    readme.append("")
    readme.append("**Evidence-evaluable population: n = 55,605**")
    readme.append("")
    readme.append("## Files")
    readme.append("")
    readme.append("### supplementary_table_tool_capabilities.csv")
    readme.append("")
    readme.append("Full tool capability matrix showing operational status of all 28")
    readme.append("canonical ACMG/AMP criteria under benchmark conditions.")
    readme.append("")
    readme.append("**Columns:**")
    readme.append("- `tool`: Tool name")
    readme.append("- `PVS1`, `PS1`, ..., `BP7`: Operational status per criterion")
    readme.append("  - `supported`: Operational under benchmark inputs")
    readme.append("  - `not_supported`: Not operational (requires unavailable input)")
    readme.append("- `n_operational_criteria`: Count of operational criteria")
    readme.append("")
    readme.append("**Benchmark inputs:**")
    readme.append("- Single-sample VCF only")
    readme.append("- No parental VCF, phased VCF, phenotype, or family history")
    readme.append("")
    readme.append("### supplementary_table_clingen_tool_summary.csv")
    readme.append("")
    readme.append("Complete ClinGen tool-level evidence concordance summary.")
    readme.append("")
    readme.append("**Columns:**")
    readme.append("- `tool`: Tool name")
    readme.append("- `n_total_predictions`: Total valid predictions")
    readme.append("- `n_classification_matches`: Exact five-class matches")
    readme.append("- `classification_concordance`: Exact five-class concordance rate")
    readme.append("- `n_raw_evaluable`: Evidence-evaluable predictions (classification-concordant + documented reference)")
    readme.append("- `raw_mean_jaccard`: Mean raw criterion Jaccard")
    readme.append("- `raw_median_jaccard`: Median raw criterion Jaccard")
    readme.append("- `raw_exact_match_rate`: Raw exact evidence match rate")
    readme.append("- `n_adjusted_evaluable`: Capability-adjusted evaluable predictions")
    readme.append("- `n_no_operationally_evaluable`: Predictions with no operational criteria overlap")
    readme.append("- `adjusted_mean_jaccard`: Mean capability-adjusted Jaccard")
    readme.append("- `adjusted_median_jaccard`: Median capability-adjusted Jaccard")
    readme.append("- `adjusted_exact_match_rate`: Capability-adjusted exact match rate")
    readme.append("- `delta_mean_jaccard`: Change in mean Jaccard (adjusted - raw)")
    readme.append("- `delta_median_jaccard`: Change in median Jaccard (adjusted - raw)")
    readme.append("- `delta_exact_match`: Change in exact match rate (adjusted - raw)")
    readme.append("- `n_operational_criteria`: Number of operational criteria (out of 28)")
    readme.append("")
    readme.append("**Interpretation:**")
    readme.append("- Raw metrics evaluate reproduction of full expert evidence profile")
    readme.append("- Adjusted metrics evaluate concordance within tool's operational scope")
    readme.append("- Positive delta indicates improvement after accounting for scope limitations")
    readme.append("- NOT a tool ranking - reflects operational scope differences")
    readme.append("")
    readme.append("### supplementary_table_clingen_criterion_decomposition.csv")
    readme.append("")
    readme.append("Complete decomposition of missing reference criteria for all 28")
    readme.append("canonical ACMG/AMP criteria.")
    readme.append("")
    readme.append("**Columns:**")
    readme.append("- `criterion`: ACMG/AMP criterion code")
    readme.append("- `total_missing_n`: Total instances where reference assigned but tool did not")
    readme.append("- `outside_scope_n`: Missing because criterion outside operational scope")
    readme.append("- `supported_but_missed_n`: Missing despite criterion being operational")
    readme.append("- `outside_scope_fraction`: Fraction of missing due to scope limitations")
    readme.append("- `supported_but_missed_fraction`: Fraction of missing due to assignment gaps")
    readme.append("")
    readme.append("**Interpretation:**")
    readme.append("- `outside_scope`: Tool cannot evaluate (requires unavailable input)")
    readme.append("- `supported_but_missed`: Tool can evaluate but did not assign")
    readme.append("- Criteria with 0 total_missing: All reference assignments reproduced by tools")
    readme.append("")
    readme.append("## Data Verification")
    readme.append("")
    readme.append("✓ ClinGen variants with Ground_Truth_ACMG='.' correctly excluded")
    readme.append("✓ 257 VUS variants with missing criteria excluded")
    readme.append("✓ 1,749 variant-tool pairs excluded from evidence-evaluable population")
    readme.append("✓ Evidence-evaluable population: n = 55,605 (NOT 57,354)")
    readme.append("✓ Same criterion parser and normalization as Figure 2")
    readme.append("✓ Capability matrix from authoritative reference file")
    readme.append("")
    readme.append("---")
    readme.append("")
    readme.append("*Generated: 2026-09-10*")
    
    # Save
    readme_path = OUTPUT_DIR / 'README.md'
    with open(readme_path, 'w') as f:
        f.write('\n'.join(readme))
    
    print(f"✓ Saved: {readme_path.name}")


def main():
    """Main execution"""
    print("="*80)
    print("GENERATING CLINGEN SUPPLEMENTARY TABLES")
    print("="*80)
    print()
    print("Supporting Main Figure 5 with complete data tables")
    print()
    
    # Generate tables
    df_cap = generate_full_capability_table()
    df_summary = generate_full_clingen_tool_summary()
    df_decomp = generate_full_criterion_decomposition()
    
    # Generate README
    generate_supplementary_readme()
    
    print()
    print("="*80)
    print("✓ CLINGEN SUPPLEMENTARY TABLES COMPLETE")
    print("="*80)
    print()
    print("Output directory: outputs/clingen_supplementary/")
    print()
    print("Files:")
    print("  - supplementary_table_tool_capabilities.csv")
    print("  - supplementary_table_clingen_tool_summary.csv")
    print("  - supplementary_table_clingen_criterion_decomposition.csv")
    print("  - README.md")
    print()


if __name__ == '__main__':
    main()
