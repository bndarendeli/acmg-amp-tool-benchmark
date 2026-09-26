#!/usr/bin/env python3
"""
Generate Supplementary Figure S3 - FOXL2 Evidence Concordance

Two-panel dumbbell plot showing raw vs capability-adjusted Mean and Total
Criterion Jaccard for FOXL2 dataset.

Matches the styling of main Figure 5 panels d/e.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# Set publication-quality style
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# Semantic colors (matching main Figure 5)
SEMANTIC_COLORS = {
    'raw': '#3d3d3d',           # Dark charcoal for raw metrics
    'adjusted': '#4A9B8E',      # Muted teal for adjusted/within-scope
    'connector': '#d8d8d8',     # Very light gray for connector lines
}

# Tool display names
TOOL_NAMES = {
    'InterVar_2018': 'InterVar 2018',
    'InterVar_2025': 'InterVar 2025',
    'BIAS': 'BIAS',
    'DiabloACMG': 'DiabloACMG',
    'Genebe': 'Genebe',
    'TAPES': 'TAPES',
    'Franklin': 'Franklin',
    'AutoGVP': 'AutoGVP',
    'VIP-HL': 'VIP-HL',
}


def load_data():
    """Load FOXL2 Total Criterion Jaccard data"""
    base_dir = Path(__file__).parent
    output_dir = base_dir / 'outputs' / 'foxl2_supplementary'
    
    # Load Total CJ data (includes both Mean and Total CJ)
    df = pd.read_csv(output_dir / 'foxl2_total_criterion_jaccard.csv')
    
    return df


def validate_data(df):
    """Validate data before plotting"""
    print("\nValidating data...")
    
    errors = []
    
    # Check number of tools
    n_tools = len(df)
    print(f"  Number of tools: {n_tools}")
    if n_tools != 9:
        errors.append(f"Expected 9 FOXL2 tools, found {n_tools}")
    
    # Check all required columns present
    required_cols = [
        'tool', 'raw_mean_jaccard', 'adjusted_mean_jaccard',
        'raw_total_criterion_jaccard', 'adjusted_total_criterion_jaccard'
    ]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        errors.append(f"Missing columns: {missing_cols}")
    
    # Check value ranges
    for col in ['raw_mean_jaccard', 'adjusted_mean_jaccard', 
                'raw_total_criterion_jaccard', 'adjusted_total_criterion_jaccard']:
        if col in df.columns:
            if df[col].min() < 0 or df[col].max() > 1:
                errors.append(f"{col} contains values outside [0, 1]")
    
    if errors:
        print("  ✗ Validation errors:")
        for error in errors:
            print(f"    - {error}")
        raise ValueError("Data validation failed")
    else:
        print("  ✓ All validation checks passed")
    
    return True


def prepare_data(df):
    """Prepare data for plotting with consistent tool ordering"""
    print("\nPreparing data...")
    
    # Add tool display names
    df['tool_name'] = df['tool'].map(TOOL_NAMES)
    
    # Sort by raw Mean CJ (descending) - highest at top
    df = df.sort_values('raw_mean_jaccard', ascending=False).reset_index(drop=True)
    
    print(f"  Tool ordering (by raw Mean CJ, descending):")
    for i, row in df.iterrows():
        print(f"    {i+1}. {row['tool_name']}: {row['raw_mean_jaccard']:.4f}")
    
    return df


def plot_dumbbell_panel(ax, df, raw_col, adj_col, x_label, panel_label, show_legend=False):
    """
    Plot dumbbell panel with semantic colors.
    
    Matches the styling of main Figure 5 panels d/e.
    Franklin is excluded from capability-adjusted display (raw value only).
    """
    
    n_tools = len(df)
    y_positions = np.arange(n_tools)
    
    # Separate Franklin from other tools
    franklin_mask = df['tool_name'] == 'Franklin'
    non_franklin_mask = ~franklin_mask
    
    # Plot connector lines (only for non-Franklin tools)
    for i, (_, row) in enumerate(df.iterrows()):
        if row['tool_name'] != 'Franklin':
            raw_val = row[raw_col]
            adj_val = row[adj_col]
            ax.plot([raw_val, adj_val], [i, i],
                   color=SEMANTIC_COLORS['connector'], linewidth=1.5, zorder=1)
    
    # Plot raw points (all tools including Franklin)
    ax.scatter(df[raw_col], y_positions,
              color=SEMANTIC_COLORS['raw'], s=50, alpha=0.9,
              edgecolor='white', linewidth=0.5, zorder=3, label='Raw')
    
    # Plot adjusted points (only non-Franklin tools)
    ax.scatter(df.loc[non_franklin_mask, adj_col], y_positions[non_franklin_mask],
              color=SEMANTIC_COLORS['adjusted'], s=50, alpha=0.9,
              edgecolor='white', linewidth=0.5, zorder=3, label='Capability-adjusted')
    
    # Y-axis labels (tool names)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(df['tool_name'], fontsize=7)
    
    # X-axis
    ax.set_xlabel(x_label, fontsize=8, fontweight='normal')
    ax.set_xlim(-0.05, 1.05)
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', length=0)
    
    # Grid
    ax.grid(True, axis='x', alpha=0.2, linewidth=0.3)
    ax.set_axisbelow(True)
    
    # Y-axis line visible (matching main Figure 5)
    ax.spines['left'].set_visible(True)
    
    # Panel label
    ax.text(-0.12, 1.05, panel_label, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left')
    
    # Legend (only if requested)
    if show_legend:
        ax.legend(loc='upper right', fontsize=6, frameon=True, 
                 edgecolor='gray', framealpha=0.9)


def generate_figure(df, output_dir):
    """Generate the two-panel Supplementary Figure S3"""
    
    print("\n" + "="*80)
    print("GENERATING SUPPLEMENTARY FIGURE S3")
    print("="*80)
    
    # Create figure with two panels
    fig = plt.figure(figsize=(8, 5))
    
    # GridSpec: 1 row × 5 columns
    # Panel a: cols 0-1, whitespace: col 2, Panel b: cols 3-4
    gs = GridSpec(1, 5, figure=fig,
                  width_ratios=[2.0, 0.8, 2.0, 0.8, 0.2],
                  wspace=0.05,
                  left=0.12, right=0.98, top=0.95, bottom=0.12)
    
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 2])
    
    print("\nPlotting panels...")
    
    # Panel a: Raw vs Adjusted Mean CJ (with legend)
    print("  ✓ Panel a: Raw vs Adjusted Mean Criterion Jaccard")
    plot_dumbbell_panel(ax_a, df,
                       'raw_mean_jaccard', 'adjusted_mean_jaccard',
                       'Mean criterion Jaccard',
                       'a', show_legend=True)
    
    # Panel b: Raw vs Adjusted Total CJ (no legend - shared with a)
    print("  ✓ Panel b: Raw vs Adjusted Total Criterion Jaccard")
    plot_dumbbell_panel(ax_b, df,
                       'raw_total_criterion_jaccard', 'adjusted_total_criterion_jaccard',
                       'Total criterion Jaccard',
                       'b', show_legend=False)
    
    # Save figure
    print("\nSaving figure...")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    pdf_path = output_dir / 'Supplementary_Figure_S3.pdf'
    png_path = output_dir / 'Supplementary_Figure_S3.png'
    
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    print(f"  ✓ Saved: {pdf_path}")
    print(f"  ✓ Saved: {png_path}")
    
    plt.close()
    
    return df


def generate_report(df):
    """Generate final report"""
    
    print("\n" + "="*80)
    print("SUPPLEMENTARY FIGURE S3 GENERATION REPORT")
    print("="*80)
    
    print("\n1. TOOL ORDERING (by raw Mean CJ, descending):")
    for i, row in df.iterrows():
        print(f"   {i+1}. {row['tool_name']}")
    
    print("\n2. OUTPUT PATHS:")
    print("   ✓ outputs/foxl2_supplementary/Supplementary_Figure_S3.pdf")
    print("   ✓ outputs/foxl2_supplementary/Supplementary_Figure_S3.png")
    
    print("\n3. DATA VERIFICATION:")
    print(f"   ✓ Number of tools: {len(df)}")
    print(f"   ✓ All tools present in both panels")
    print(f"   ✓ Identical tool ordering in both panels")
    print(f"   ✓ All values in range [0, 1]")
    
    print("\n4. METRICS SUMMARY:")
    print(f"\n   {'Tool':<20} {'Raw Mean':>10} {'Adj Mean':>10} {'Raw Total':>10} {'Adj Total':>10}")
    print("   " + "-"*70)
    for _, row in df.iterrows():
        print(f"   {row['tool_name']:<20} {row['raw_mean_jaccard']:>10.4f} {row['adjusted_mean_jaccard']:>10.4f} {row['raw_total_criterion_jaccard']:>10.4f} {row['adjusted_total_criterion_jaccard']:>10.4f}")
    
    print("\n" + "="*80)
    print("✓ SUPPLEMENTARY FIGURE S3 GENERATION COMPLETE")
    print("="*80)


def main():
    """Main execution"""
    print("="*80)
    print("SUPPLEMENTARY FIGURE S3 - FOXL2 EVIDENCE CONCORDANCE")
    print("="*80)
    print()
    print("Two-panel dumbbell plot:")
    print("  Panel a: Raw vs Capability-adjusted Mean Criterion Jaccard")
    print("  Panel b: Raw vs Capability-adjusted Total Criterion Jaccard")
    print()
    
    # Load data
    print("Loading data...")
    df = load_data()
    print(f"  ✓ Loaded {len(df)} tools")
    
    # Validate data
    validate_data(df)
    
    # Prepare data with consistent ordering
    df = prepare_data(df)
    
    # Generate figure
    output_dir = Path(__file__).parent / 'outputs' / 'foxl2_supplementary'
    df = generate_figure(df, output_dir)
    
    # Generate report
    generate_report(df)


if __name__ == '__main__':
    main()
