#!/usr/bin/env python3
"""
Generate Refined Figure 5 - Publication Quality

Three-panel design with restrained Nature Communications styling:
- Panel A: Tool-level scatter with tool-specific colors and direct labeling
- Panel B: Raw vs adjusted dumbbell with semantic colors
- Panel C: Mechanisms with semantic colors matching Panel B

NO panel titles - meaning conveyed through axes and visual encoding.
Evidence-evaluable population: n = 55,605
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from pathlib import Path
from adjustText import adjust_text
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

# Semantic colors for Panels B and C
SEMANTIC_COLORS = {
    'raw': '#3d3d3d',           # Dark charcoal for raw metrics
    'adjusted': '#4A9B8E',      # Muted teal for adjusted/within-scope
    'outside_scope': '#a8a8a8', # Medium gray for outside scope
    'connector': '#d8d8d8',     # Very light gray for connector lines
}

# Tool-specific colors for Panel A only
# Muted, colorblind-friendly qualitative palette
TOOL_COLORS = {
    'Genebe': '#A23B72',        # Blue
    'Franklin': '#DE8F05',      # Orange
    'CPSR': '#029E73',          # Green
    'VIP-HL': '#CC78BC',        # Purple
    'CancerSIGVAR': '#CA9161',  # Tan
    'InterVar_2018': '#949494', # Gray
    'BIAS': '#ECE133',          # Yellow
    'CharGer_Local': '#56B4E9', # Light blue
    'CharGer_Online': '#009E73', # Teal
    'TAPES': '#F0E442',         # Light yellow
    'AutoGVP': '#0072B2',       # Dark blue
    'DiabloACMG': '#D55E00',    # Red-orange
    'Exomiser': '#CC79A7',      # Pink
    'InterVar_2025': '#7B61A8', # Light gray
}

TOOL_NAMES = {
    'InterVar_2018': 'InterVar 2018',
    'InterVar_2025': 'InterVar 2025',
    'BIAS': 'BIAS',
    'CharGer_Local': 'CharGer Local',
    'CharGer_Online': 'CharGer Online',
    'DiabloACMG': 'DiabloACMG',
    'Exomiser': 'Exomiser',
    'Genebe': 'Genebe',
    'TAPES': 'TAPES',
    'Franklin': 'Franklin',
    'AutoGVP': 'AutoGVP',
    'VIP-HL': 'VIP-HL',
    'CancerSIGVAR': 'CancerSIGVAR',
    'CPSR': 'CPSR'
}

def load_data():
    """Load all required data files"""
    data = {}
    data['evidence'] = pd.read_csv('outputs/fig5_evidence_distributions.csv')
    data['master'] = pd.read_csv('outputs/fig5_master_variant_tool.csv')
    data['adjusted'] = pd.read_csv('outputs/fig5_capability_adjusted_summary.csv')
    data['capability'] = pd.read_csv('outputs/tool_criterion_capability_matrix.csv')
    data['decomposition'] = pd.read_csv('outputs/fig5_missing_criteria_decomposition.csv')
    return data


def prepare_panel_a_data(data):
    """Prepare data for Panel A"""
    df_evidence = data['evidence']
    df_master = data['master']
    
    df_clingen_ev = df_evidence[
        (df_evidence['dataset'] == 'clingen_28012026') &
        (df_evidence['tool'] != 'POOLED')
    ].copy()
    
    df_clingen_master = df_master[df_master['dataset'] == 'clingen_28012026'].copy()
    
    results = []
    for tool in df_clingen_ev['tool'].unique():
        df_tool = df_clingen_master[df_clingen_master['tool'] == tool]
        n_total = len(df_tool)
        n_concordant = len(df_tool[df_tool['classification_match'] == True])
        concordance = n_concordant / n_total if n_total > 0 else 0
        
        ev_row = df_clingen_ev[df_clingen_ev['tool'] == tool].iloc[0]
        median_j = ev_row['median_criterion_jaccard']
        
        results.append({
            'tool': tool,
            'tool_name': TOOL_NAMES.get(tool, tool),
            'classification_concordance': concordance,
            'median_jaccard': median_j,
        })
    
    return pd.DataFrame(results)


def prepare_panel_b_data(data):
    """Prepare data for Panel B with operational coverage"""
    df_adjusted = data['adjusted'].copy()
    df_capability = data['capability'].copy()
    
    coverage = []
    for tool in df_capability['tool'].unique():
        df_tool = df_capability[df_capability['tool'] == tool]
        n_operational = len(df_tool[df_tool['operational_status'] == 'supported'])
        coverage.append({'tool': tool, 'n_operational': n_operational})
    
    df_coverage = pd.DataFrame(coverage)
    df_panel_b = df_adjusted.merge(df_coverage, on='tool')
    df_panel_b['tool_name'] = df_panel_b['tool'].map(TOOL_NAMES)
    
    # Sort by raw median Jaccard (descending)
    df_panel_b = df_panel_b.sort_values('raw_median_jaccard', ascending=False)
    
    return df_panel_b


def prepare_panel_c_data(data):
    """Prepare data for Panel C"""
    df_decomp = data['decomposition'].copy()
    df_pooled = df_decomp[df_decomp['analysis_type'] == 'pooled'].copy()
    
    criterion_order = ['PM3', 'PP1', 'PP4', 'PS4', 'BP4']
    df_panel_c = df_pooled[df_pooled['criterion'].isin(criterion_order)].copy()
    
    df_panel_c['criterion'] = pd.Categorical(df_panel_c['criterion'], 
                                             categories=criterion_order, 
                                             ordered=True)
    df_panel_c = df_panel_c.sort_values('criterion')
    
    return df_panel_c


def plot_panel_a(ax, df):
    """Plot Panel A: highlighted tools + contextual gray points"""

    # Tools selected to illustrate informative patterns
    HIGHLIGHT_TOOLS = {
        'Genebe',
        'AutoGVP',
        'Franklin',
        'TAPES',
        'VIP-HL',
        'CPSR',
        'InterVar_2025'
    }

    # Manual label offsets for highlighted tools only
    LABEL_OFFSETS = {
        'Genebe':         (-8, -6),
        'AutoGVP':        (-10, -7),
        'Franklin':       (0, 9),
        'TAPES':  (5, -9),
        'VIP-HL':         (-8, -3),
        'CPSR':           (7, 4),
        'InterVar_2025':  (4, -7),
    }

    # ---------------------------------------------------------
    # 1. Plot all tools
    # ---------------------------------------------------------
    for _, row in df.iterrows():

        if row['tool'] in HIGHLIGHT_TOOLS:
            # Highlighted tools retain their tool-specific color
            color = TOOL_COLORS.get(row['tool'], '#666666')
            alpha = 0.95
            point_size = 48
            zorder = 4
        else:
            # Remaining tools provide context only
            color = '#b8b8b8'
            alpha = 0.75
            point_size = 36
            zorder = 2

        ax.scatter(
            row['classification_concordance'],
            row['median_jaccard'],
            s=point_size,
            color=color,
            alpha=alpha,
            edgecolors='white',
            linewidth=0.5,
            zorder=zorder
        )

    # ---------------------------------------------------------
    # 2. Label ONLY highlighted tools
    # ---------------------------------------------------------
    highlighted_df = df[df['tool'].isin(HIGHLIGHT_TOOLS)]

    for _, row in highlighted_df.iterrows():

        dx, dy = LABEL_OFFSETS.get(row['tool'], (5, 5))

        ax.annotate(
            row['tool_name'],
            xy=(
                row['classification_concordance'],
                row['median_jaccard']
            ),
            xytext=(dx, dy),
            textcoords='offset points',
            fontsize=6.3,
            color='#333333',
            ha='left' if dx >= 0 else 'right',
            va='center',
            zorder=5
        )

    # ---------------------------------------------------------
    # 3. Axes
    # ---------------------------------------------------------
    ax.set_xlabel(
        'Classification concordance (valid predictions)',
        fontsize=8
    )

    ax.set_ylabel(
        'Median criterion Jaccard',
        fontsize=8
    )

    ax.set_xlim(0.30, 1.02)
    ax.set_ylim(0.15, 1.05)

    # Very subtle grid
    ax.grid(
        True,
        alpha=0.08,
        linestyle='-',
        linewidth=0.3,
        color='#d8d8d8',
        zorder=0
    )

    ax.set_axisbelow(True)

    # Clean publication-style spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)

    # Panel label
    ax.text(
        -0.20,
        1.08,
        'a',
        transform=ax.transAxes,
        fontsize=10,
        fontweight='bold',
        va='top',
        ha='left'
    )


def plot_panel_b(ax, df):
    """Plot Panel B: Raw vs adjusted with semantic colors"""

    n_tools = len(df)
    y_positions = np.arange(n_tools)

    # Plot dumbbells
    for i, (_, row) in enumerate(df.iterrows()):

        raw = row['raw_median_jaccard']
        adjusted = row['adjusted_median_jaccard']

        # Connector
        ax.plot(
            [raw, adjusted],
            [i, i],
            color=SEMANTIC_COLORS['connector'],
            linewidth=1.0,
            alpha=0.7,
            zorder=1
        )

        # Same raw and adjusted value
        if np.isclose(raw, adjusted):

            # Raw = larger outer charcoal circle
            ax.scatter(
                raw, i,
                s=60,
                color=SEMANTIC_COLORS['raw'],
                edgecolors='white',
                linewidth=0.6,
                zorder=3
            )

            # Adjusted = smaller teal circle on top
            ax.scatter(
                adjusted, i,
                s=30,
                color=SEMANTIC_COLORS['adjusted'],
                edgecolors='white',
                linewidth=0.4,
                zorder=4
            )

        else:

            # Raw
            ax.scatter(
                raw, i,
                s=50,
                color=SEMANTIC_COLORS['raw'],
                edgecolors='white',
                linewidth=0.6,
                zorder=3
            )

            # Adjusted
            ax.scatter(
                adjusted, i,
                s=50,
                color=SEMANTIC_COLORS['adjusted'],
                edgecolors='white',
                linewidth=0.6,
                zorder=3
            )

    # Y-axis labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels(df['tool_name'], fontsize=7)
    ax.invert_yaxis()

    # Operational criteria header
    ax.text(
        1.10,
        1.02,
        'Operational\ncriteria',
        va='bottom',
        ha='center',
        fontsize=6.5,
        color='#666666',
        linespacing=0.9,
        transform=ax.transAxes,
        clip_on=False
    )

    # Operational coverage + small n
    for i, (_, row) in enumerate(df.iterrows()):

        coverage_text = f"{row['n_operational']}/28"

        ax.text(
            1.10,
            i,
            coverage_text,
            va='center',
            ha='center',
            fontsize=6.5,
            color='#666666',
            transform=ax.get_yaxis_transform()
        )

        if row['n_evaluable'] < 100:
            ax.text(
                1.16,
                i,
                f"n={int(row['n_evaluable'])}",
                va='center',
                ha='left',
                fontsize=5.5,
                color='#999999',
                style='italic',
                transform=ax.get_yaxis_transform()
            )

    # X-axis
    ax.set_xlabel('Median criterion Jaccard', fontsize=8)
    ax.set_xlim(-0.05, 1.05)

    # Legend
    legend_elements = [
        plt.Line2D(
            [0], [0],
            marker='o',
            color='w',
            markerfacecolor=SEMANTIC_COLORS['raw'],
            markersize=6,
            markeredgecolor='white',
            markeredgewidth=0.6,
            label='Raw',
            linestyle='None'
        ),
        plt.Line2D(
            [0], [0],
            marker='o',
            color='w',
            markerfacecolor=SEMANTIC_COLORS['adjusted'],
            markersize=6,
            markeredgecolor='white',
            markeredgewidth=0.6,
            label='Adjusted',
            linestyle='None'
        )
    ]

    ax.legend(
        handles=legend_elements,
        loc='lower center',
        bbox_to_anchor=(0.50, 1.01),
        ncol=2,
        fontsize=6.5,
        frameon=False,
        handletextpad=0.3,
        columnspacing=1.2
    )

    # Minimal gridlines
    ax.grid(
        True,
        axis='x',
        alpha=0.12,
        linestyle='-',
        linewidth=0.3,
        color='#d8d8d8',
        zorder=0
    )

    ax.set_axisbelow(True)

    # Clean spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)

    # Panel label
    ax.text(
        -0.14,
        1.08,
        'b',
        transform=ax.transAxes,
        fontsize=10,
        fontweight='bold',
        va='top',
        ha='left'
    )

def plot_panel_c(ax, df):
    """Plot Panel C: Mechanisms with semantic colors"""
    
    n_criteria = len(df)
    y_positions = np.arange(n_criteria)
    
    # Plot stacked bars with semantic colors
    outside_pct = df['outside_scope_fraction_of_missing'].values * 100
    ax.barh(y_positions, outside_pct, height=0.60,
           color=SEMANTIC_COLORS['outside_scope'], edgecolor='white', 
           linewidth=0.5)
    
    supported_pct = df['supported_but_missed_fraction_of_missing'].values * 100
    ax.barh(y_positions, supported_pct, left=outside_pct,
           height=0.60, color=SEMANTIC_COLORS['adjusted'], 
           edgecolor='white', linewidth=0.5)
    
    # Add percentage labels inside bars
    for i, (_, row) in enumerate(df.iterrows()):
        outside = row['outside_scope_fraction_of_missing'] * 100
        supported = row['supported_but_missed_fraction_of_missing'] * 100
        
        # Outside scope label (if > 8%)
        if outside > 8:
            ax.text(outside / 2, i, f"{outside:.0f}%",
                   va='center', ha='center', fontsize=6.5, 
                   color='white', fontweight='normal')
        
        # Supported but missed label (if > 8%)
        if supported > 8:
            ax.text(outside + supported / 2, i, f"{supported:.0f}%",
                   va='center', ha='center', fontsize=6.5, 
                   color='white', fontweight='normal')
        elif supported > 3:  # Show small PP4 segment
            ax.text(outside + supported / 2, i, f"{supported:.0f}%",
                   va='center', ha='center', fontsize=5.5, 
                   color='white', fontweight='normal')
    
    # Y-axis labels (criterion names only)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(df['criterion'].values, fontsize=7)
    ax.invert_yaxis()
    
    # Add total counts on the right
    for i, (_, row) in enumerate(df.iterrows()):
        count_text = f"n = {row['n_total_missing']:,}"
        ax.text(104, i, count_text, va='center', ha='left',
               fontsize=6, color='#666666')
    
    # X-axis
    ax.set_xlabel('Missing reference criterion assignments (%)', 
                  fontsize=8)
    ax.set_xlim(0, 118)
    
    # Legend - place to the right
    legend_elements = [
        mpatches.Patch(facecolor=SEMANTIC_COLORS['outside_scope'], 
                      edgecolor='white', linewidth=0.5,
                      label='Outside operational scope'),
        mpatches.Patch(facecolor=SEMANTIC_COLORS['adjusted'], 
                      edgecolor='white', linewidth=0.5,
                      label='Operationally evaluable but not assigned')
    ]
    ax.legend(handles=legend_elements,
    loc='upper center',
    bbox_to_anchor=(0.5, -0.23),
    ncol=2,
    fontsize=6.5,
    frameon=False,
    handletextpad=0.5,
    columnspacing=1.5)
    
    # Minimal gridlines
    ax.grid(True, axis='x', alpha=0.12, linestyle='-', 
           linewidth=0.3, color='#d8d8d8', zorder=0)
    ax.set_axisbelow(True)
    
    # Clean spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    
    # Panel label
    ax.text(-0.09, 1.08, 'c', transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left')


def generate_refined_figure(data, output_dir):
    """Generate refined publication-quality Figure 5"""
    print("\n" + "="*80)
    print("GENERATING REFINED FIGURE 4 - PUBLICATION QUALITY")
    print("="*80)
    
    # Prepare data
    df_a = prepare_panel_a_data(data)
    df_b = prepare_panel_b_data(data)
    df_c = prepare_panel_c_data(data)
    
    # Create figure with refined layout
    fig = plt.figure(figsize=(7.1, 7.0))
    
    # Create custom grid with reduced vertical spacing
    gs = GridSpec(2, 2, figure=fig, 
                  height_ratios=[1, 0.42],
                  width_ratios=[0.43, 0.57],
                  hspace=0.32, wspace=0.45,
                  left=0.08, right=0.92, top=0.97, bottom=0.12)
    
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    
    # Plot panels
    print("\nPlotting panels...")
    plot_panel_a(ax_a, df_a)
    print("✓ Panel A plotted")
    plot_panel_b(ax_b, df_b)
    print("✓ Panel B plotted")
    plot_panel_c(ax_c, df_c)
    print("✓ Panel C plotted")
    
    # Save preview
    print("\n" + "="*80)
    print("SAVING PREVIEW")
    print("="*80)
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    preview_path = output_dir / 'Fig5_Evidence_Concordance.png'
    fig.savefig(preview_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Preview saved: {preview_path.name}")
    pdf_path = output_dir / "Fig5_Evidence_Concordance.pdf"
    fig.savefig(pdf_path, format="pdf", dpi=300, bbox_inches="tight", pad_inches=0.02)
    print(f"✓ Preview saved: {pdf_path.name}")

    plt.close(fig)
    
    return df_a, df_b, df_c


def main():
    """Main execution"""
    print("="*80)
    print("FIGURE 5 REFINEMENT - PUBLICATION QUALITY")
    print("="*80)
    print()
    print("Refinements:")
    print("  - No panel titles")
    print("  - Panel A: Tool-specific colors, all tools labeled")
    print("  - Panels B & C: Semantic colors (charcoal/teal)")
    print("  - Minimal typography")
    print("  - Reduced spacing")
    print()
    
    data = load_data()
    output_dir = Path('outputs/figures')
    df_a, df_b, df_c = generate_refined_figure(data, output_dir)
    #generate_refinement_report(df_a, df_b, df_c, output_dir)
    
    print()
    print("="*80)
    print("✓ REFINEMENT COMPLETE")
    print("="*80)
    print()
    print("Preview file:")
    print("  - outputs/figures/fig5_redesign_preview_v2.png")
    print()
    print("Review this preview before generating final SVG/PDF versions.")
    print()


if __name__ == '__main__':
    main()
