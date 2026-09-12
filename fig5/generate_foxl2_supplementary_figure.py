#!/usr/bin/env python3
"""
Generate FOXL2 Supplementary Figure: Evidence Concordance Analysis

Creates a 2-panel horizontal figure showing:
- Panel a: Raw vs capability-adjusted evidence concordance (dumbbell plot)
- Panel b: Decomposition of missing reference evidence (stacked bars)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Configuration
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'outputs' / 'figures'
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Input files
SUMMARY_FILE = BASE_DIR / 'outputs' / 'foxl2_supplementary' / 'foxl2_capability_adjusted_summary.csv'
DECOMP_FILE = BASE_DIR / 'outputs' / 'foxl2_supplementary' / 'foxl2_criterion_missing_decomposition.csv'

# Output files
OUTPUT_PDF = OUTPUT_DIR / 'FOXL2_Supplementary_Evidence_Concordance.pdf'
OUTPUT_PNG = OUTPUT_DIR / 'FOXL2_Supplementary_Evidence_Concordance.png'

# Display name mapping
TOOL_DISPLAY_NAMES = {
    'Genebe': 'GeneBe',
    'InterVar_2018': 'InterVar 2018',
    'InterVar_2025': 'InterVar 2025',
    'CharGer_Local': 'CharGer (local)',
    'CharGer_Online': 'CharGer (online)',
    'BIAS': 'BIAS',
    'DiabloACMG': 'DiabloACMG',
    'Franklin': 'Franklin',
    'AutoGVP': 'AutoGVP',
    'TAPES': 'TAPES',
    'VIP-HL': 'VIP-HL',
    'CPSR': 'CPSR',
    'Exomiser': 'Exomiser',
    'CancerSIGVAR': 'CancerSIGVAR'
}

# Colors (consistent with Main Figure 5)
COLORS = {
    'raw': '#2D2D2D',              # Black for raw
    'adjusted': '#4A9B9B',          # Teal for adjusted
    'outside_scope': '#999999',     # Grey for outside scope
    'connector': '#CCCCCC'          # Light grey for connector lines
}


def load_data():
    """Load summary and decomposition data"""
    print("="*80)
    print("LOADING DATA")
    print("="*80)
    
    df_summary = pd.read_csv(SUMMARY_FILE)
    df_decomp = pd.read_csv(DECOMP_FILE)
    
    print(f"\nSummary data: {len(df_summary)} tools")
    print(f"Decomposition data: {len(df_decomp)} criteria")
    
    return df_summary, df_decomp


def create_panel_a(ax, df_summary):
    """Create Panel a: Raw vs capability-adjusted evidence concordance"""
    print("\n" + "="*80)
    print("PANEL A: RAW VS CAPABILITY-ADJUSTED CONCORDANCE")
    print("="*80)
    
    # Prepare data
    df = df_summary.copy()
    
    # Apply display names
    df['tool_display'] = df['tool'].map(TOOL_DISPLAY_NAMES)
    
    # Sort by decreasing capability-adjusted median Jaccard
    df = df.sort_values('adjusted_median_jaccard', ascending=True)
    
    # Print values
    print("\nTool-level concordance:")
    print(f"{'Tool':<20} {'Raw Median':>12} {'Adj Median':>12} {'Delta':>10} {'n':>6}")
    print("-" * 70)
    for _, row in df.iterrows():
        delta = row['adjusted_median_jaccard'] - row['raw_median_jaccard']
        print(f"{row['tool_display']:<20} {row['raw_median_jaccard']:>12.3f} "
              f"{row['adjusted_median_jaccard']:>12.3f} {delta:>+10.3f} {row['n_adjusted_evaluable']:>6}")
    
    # Plot
    n_tools = len(df)
    y_pos = np.arange(n_tools)
    
    # Connector lines
    for i, (_, row) in enumerate(df.iterrows()):
        ax.plot([row['raw_median_jaccard'], row['adjusted_median_jaccard']], 
                [i, i], 
                color=COLORS['connector'], 
                linewidth=1.5, 
                zorder=1)
    
    # Raw values (black circles)
    ax.scatter(df['raw_median_jaccard'], y_pos,
              s=80, c=COLORS['raw'], marker='o',
              edgecolors='white', linewidths=0.5,
              zorder=3, label='Raw')
    
    # Adjusted values (teal circles)
    ax.scatter(df['adjusted_median_jaccard'], y_pos,
              s=80, c=COLORS['adjusted'], marker='o',
              edgecolors='white', linewidths=0.5,
              zorder=3, label='Capability-adjusted')
    
    # Add n values to the right
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(1.02, i, f"n = {row['n_adjusted_evaluable']}", 
                va='center', ha='left', fontsize=8, color='#666666')
    
    # Styling
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['tool_display'], fontsize=9)
    ax.set_xlim(-0.05, 1.15)
    ax.set_xlabel('Median criterion Jaccard similarity', fontsize=9, labelpad=5)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])
    
    # Grid
    ax.grid(axis='x', alpha=0.2, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', 
                  markerfacecolor=COLORS['raw'], markersize=8,
                  markeredgecolor='white', markeredgewidth=0.5,
                  label='Raw'),
        plt.Line2D([0], [0], marker='o', color='w',
                  markerfacecolor=COLORS['adjusted'], markersize=8,
                  markeredgecolor='white', markeredgewidth=0.5,
                  label='Capability-adjusted')
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.50, -0.14),
        ncol=2,
        frameon=False,
        fontsize=8,
        handletextpad=0.4,
        columnspacing=1.0
    )
    
    # Panel label
    ax.text(-0.15, 1.05, 'a', transform=ax.transAxes,
           fontsize=12, fontweight='bold', va='top')


def create_panel_b(ax, df_decomp):
    """Create Panel b: Decomposition of missing reference evidence"""
    print("\n" + "="*80)
    print("PANEL B: MISSING REFERENCE EVIDENCE DECOMPOSITION")
    print("="*80)
    
    # Filter to criteria with missing instances
    df = df_decomp[df_decomp['n_total_missing'] > 0].copy()
    
    # Sort by decreasing total missing
    df = df.sort_values('n_total_missing', ascending=True)
    
    # Print values
    print("\nCriterion-level decomposition:")
    print(f"{'Criterion':<10} {'Total':>8} {'Outside':>8} {'Evaluable':>10} {'Outside %':>10} {'Evaluable %':>12}")
    print("-" * 70)
    for _, row in df.iterrows():
        print(f"{row['criterion']:<10} {row['n_total_missing']:>8} "
              f"{row['n_outside_operational_scope']:>8} "
              f"{row['n_operationally_evaluable_but_not_assigned']:>10} "
              f"{row['outside_scope_fraction_of_missing']*100:>9.1f}% "
              f"{row['operationally_evaluable_but_not_assigned_fraction']*100:>11.1f}%")
    
    # Calculate percentages
    outside_pct = df['outside_scope_fraction_of_missing'] * 100
    evaluable_pct = df['operationally_evaluable_but_not_assigned_fraction'] * 100
    
    # Plot
    n_criteria = len(df)
    y_pos = np.arange(n_criteria)
    
    # Outside scope (grey)
    ax.barh(y_pos, outside_pct, height=0.7,
           color=COLORS['outside_scope'], edgecolor='white', linewidth=0.5,
           label='Outside operational scope')
    
    # Operationally evaluable but not assigned (teal)
    ax.barh(y_pos, evaluable_pct, left=outside_pct, height=0.7,
           color=COLORS['adjusted'], edgecolor='white', linewidth=0.5,
           label='Operationally evaluable but not assigned')
    
    # Add total counts to the right
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(102, i, f"n = {row['n_total_missing']}", 
                va='center', ha='left', fontsize=8, color='#666666')
    
    # Styling
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['criterion'], fontsize=9, family='monospace')
    ax.set_xlim(0, 115)
    ax.set_xlabel('Percentage of missing reference\ncriterion assignments (%)', fontsize=9, labelpad=5)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(['0', '25', '50', '75', '100'])
    
    # Grid
    ax.grid(axis='x', alpha=0.2, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['outside_scope'], 
                      edgecolor='white', linewidth=0.5,
                      label='Outside operational scope'),
        mpatches.Patch(facecolor=COLORS['adjusted'],
                      edgecolor='white', linewidth=0.5,
                      label='Operationally evaluable but not assigned')
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.50, -0.18),
        ncol=1,
        frameon=False,
        fontsize=8,
        handlelength=1.6,
        handletextpad=0.6,
        labelspacing=0.35
    )
    
    # Panel label
    ax.text(-0.15, 1.05, 'b', transform=ax.transAxes,
           fontsize=12, fontweight='bold', va='top')


def create_figure(df_summary, df_decomp):
    """Create the complete 2-panel figure"""
    print("\n" + "="*80)
    print("CREATING FIGURE")
    print("="*80)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(7.1, 4.8),
        gridspec_kw={
            'width_ratios': [1.0, 1.0],
            'wspace': 0.42
        }
    )

    # Create panels 
    create_panel_a(ax1, df_summary)
    create_panel_b(ax2, df_decomp)

    # Adjust layout
    plt.subplots_adjust(
        left=0.10,
        right=0.97,
        top=0.95,
        bottom=0.25,
        wspace=0.42
    )
    
    # Save
    print(f"\nSaving figure...")
    fig.savefig(OUTPUT_PDF, dpi=600, bbox_inches='tight')
    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches='tight')
    
    print(f"✓ Saved: {OUTPUT_PDF.name}")
    print(f"✓ Saved: {OUTPUT_PNG.name}")
    
    return fig


def print_summary_statistics(df_summary, df_decomp):
    """Print summary statistics"""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    print(f"\nPanel a (Evidence Concordance):")
    print(f"  Tools analyzed: {len(df_summary)}")
    print(f"  Total variant-tool pairs: {df_summary['n_adjusted_evaluable'].sum()}")
    print(f"  Mean raw median Jaccard: {df_summary['raw_median_jaccard'].mean():.3f}")
    print(f"  Mean adjusted median Jaccard: {df_summary['adjusted_median_jaccard'].mean():.3f}")
    print(f"  Mean delta: {(df_summary['adjusted_median_jaccard'] - df_summary['raw_median_jaccard']).mean():.3f}")
    
    df_with_missing = df_decomp[df_decomp['n_total_missing'] > 0]
    print(f"\nPanel b (Missing Evidence Decomposition):")
    print(f"  Criteria with missing instances: {len(df_with_missing)}")
    print(f"  Total missing instances: {df_with_missing['n_total_missing'].sum()}")
    print(f"  Mean outside scope fraction: {df_with_missing['outside_scope_fraction_of_missing'].mean():.3f}")
    print(f"  Mean evaluable fraction: {df_with_missing['operationally_evaluable_but_not_assigned_fraction'].mean():.3f}")
    
    # Criteria patterns
    n_100_outside = (df_with_missing['outside_scope_fraction_of_missing'] == 1.0).sum()
    n_100_evaluable = (df_with_missing['operationally_evaluable_but_not_assigned_fraction'] == 1.0).sum()
    n_mixed = len(df_with_missing) - n_100_outside - n_100_evaluable
    
    print(f"\n  100% outside scope: {n_100_outside} criteria")
    print(f"  100% evaluable but not assigned: {n_100_evaluable} criteria")
    print(f"  Mixed: {n_mixed} criteria")


def main():
    """Main execution"""
    print("="*80)
    print("FOXL2 SUPPLEMENTARY FIGURE: EVIDENCE CONCORDANCE ANALYSIS")
    print("="*80)
    print()
    print("Creating a 2-panel horizontal figure:")
    print("  Panel a: Raw vs capability-adjusted evidence concordance")
    print("  Panel b: Decomposition of missing reference evidence")
    print()
    
    # Load data
    df_summary, df_decomp = load_data()
    
    # Create figure
    fig = create_figure(df_summary, df_decomp)
    
    # Print summary
    print_summary_statistics(df_summary, df_decomp)
    
    print("\n" + "="*80)
    print("✓ FOXL2 SUPPLEMENTARY FIGURE COMPLETE")
    print("="*80)
    print()
    print(f"Output files:")
    print(f"  - {OUTPUT_PDF}")
    print(f"  - {OUTPUT_PNG}")
    print()
    print("Figure specifications:")
    print(f"  - Width: 7.1 inches (two-column)")
    print(f"  - Resolution: 600 DPI")
    print(f"  - Format: PDF (vector) + PNG (raster)")
    print(f"  - Style: Nature Communications supplementary")
    print()


if __name__ == '__main__':
    main()
