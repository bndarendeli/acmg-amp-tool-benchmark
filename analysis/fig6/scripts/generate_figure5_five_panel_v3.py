#!/usr/bin/env python3
"""
Generate Figure 5 - Five-Panel Design (v3)

Changes relative to v2:
1. Larger scatter points in panels A-C.
2. Panels D-E show the number of operationally supported ACMG/AMP
   criteria for each tool in "n/28" format.
3. Capability counts are read directly from:
       reference/tool_acmg_capabilities.csv
4. Franklin is shown as "–" because its operational criterion scope
   was not established with sufficient confidence for the
   capability-adjusted analysis.

IMPORTANT:
This is a visualization-only revision.
No performance metrics, Jaccard values, classification concordance
values, population filters, or capability-adjusted values are recalculated.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import matplotlib as mpl


# =============================================================================
# GLOBAL STYLE
# =============================================================================

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


# =============================================================================
# COLORS
# =============================================================================

SEMANTIC_COLORS = {
    'raw': '#3d3d3d',
    'adjusted': '#4A9B8E',
    'connector': '#d8d8d8',
}


TOOL_COLORS = {
    'Genebe': '#A23B72',
    'Franklin': '#DE8F05',
    'CPSR': '#029E73',
    'VIP-HL': '#CC78BC',
    'CancerSIGVAR': '#CA9161',
    'InterVar_2018': '#949494',
    'BIAS': '#ECE133',
    'CharGer_Local': '#56B4E9',
    'CharGer_Online': '#009E73',
    'TAPES': '#F0E442',
    'AutoGVP': '#0072B2',
    'DiabloACMG': '#D55E00',
    'Exomiser': '#CC79A7',
    'InterVar_2025': '#7B61A8',
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


# Highlighting in panels A-C
HIGHLIGHT_TOOLS = {
    'Genebe',
    'AutoGVP',
    'Franklin',
    'TAPES',
    'VIP-HL',
    'CPSR',
    'InterVar_2025'
}


# All standard ACMG/AMP criteria
ACMG_CRITERIA = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]

TOTAL_ACMG_CRITERIA = len(ACMG_CRITERIA)

assert TOTAL_ACMG_CRITERIA == 28


# =============================================================================
# LOAD DATA
# =============================================================================

def load_data():
    """Load all required figure data."""

    output_dir = Path('outputs')

    data = {}

    data['evidence'] = pd.read_csv(
        output_dir / 'fig4_evidence_distributions.csv'
    )

    data['master'] = pd.read_csv(
        output_dir / 'fig4_master_variant_tool.csv'
    )

    data['total_cj'] = pd.read_csv(
        output_dir / 'clingen_total_criterion_jaccard.csv'
    )

    data['adjusted'] = pd.read_csv(
        output_dir / 'fig4_capability_adjusted_summary.csv'
    )

    return data


# =============================================================================
# CAPABILITY COUNTS
# =============================================================================

def load_capability_counts():
    """
    Load the operational capability matrix that was used for the
    capability-adjusted evidence analysis.

    Source:
        outputs/tool_criterion_capability_matrix.csv

    Expected structure:
        14 tools × 28 ACMG/AMP criteria = 392 rows

    A criterion is counted only when:
        operational_status == 'supported'

    Returns:
        dict: tool -> number of operationally supported criteria
    """

    script_dir = Path(__file__).resolve().parent
    capability_path = (
        script_dir / 'outputs' / 'tool_criterion_capability_matrix.csv'
    )

    if not capability_path.exists():
        raise FileNotFoundError(
            f"Operational capability matrix not found: {capability_path}"
        )

    print("\nLoading operational capability matrix:")
    print(f"  {capability_path}")

    df = pd.read_csv(capability_path)

    print(f"  Shape: {df.shape}")
    print(f"  Columns: {list(df.columns)}")

    # --------------------------------------------------------------
    # Validate required columns
    # --------------------------------------------------------------

    required_columns = {
        'tool',
        'criterion',
        'operational_status'
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(
            "Capability matrix missing required columns: "
            + ", ".join(sorted(missing))
        )

    # --------------------------------------------------------------
    # Validate expected matrix dimensions
    # --------------------------------------------------------------

    n_tools = df['tool'].nunique()
    n_criteria = df['criterion'].nunique()

    if n_tools != 14:
        raise ValueError(
            f"Expected 14 tools, found {n_tools}"
        )

    if n_criteria != 28:
        raise ValueError(
            f"Expected 28 ACMG/AMP criteria, found {n_criteria}"
        )

    # Every tool should have exactly 28 criterion rows
    rows_per_tool = df.groupby('tool')['criterion'].nunique()

    if not (rows_per_tool == 28).all():
        raise ValueError(
            "Not every tool has exactly 28 criterion entries."
        )

    # --------------------------------------------------------------
    # Count operationally supported criteria
    # --------------------------------------------------------------

    supported = df[
        df['operational_status']
        .astype(str)
        .str.strip()
        .str.lower()
        .eq('supported')
    ].copy()

    capability_counts = (
        supported
        .groupby('tool')['criterion']
        .nunique()
        .to_dict()
    )

    # --------------------------------------------------------------
    # Validation output
    # --------------------------------------------------------------

    print("\nOperationally supported criteria:")

    for tool in sorted(df['tool'].unique()):

        n_supported = capability_counts.get(tool, 0)

        # Franklin is intentionally not capability-adjusted
        if tool == 'Franklin':
            display = '– (not adjusted)'
        else:
            display = f'{n_supported}/28'

        print(f"  {tool:<20} {display}")

    return capability_counts


# =============================================================================
# PREPARE PANELS A-C
# =============================================================================

def prepare_panel_abc_data(data):
    """Prepare data for scatter plots A-C."""

    df_evidence = data['evidence']
    df_master = data['master']
    df_total_cj = data['total_cj']

    df_clingen_ev = df_evidence[
        (df_evidence['dataset'] == 'clingen_28012026') &
        (df_evidence['tool'] != 'POOLED')
    ].copy()

    df_clingen_master = df_master[
        df_master['dataset'] == 'clingen_28012026'
    ].copy()

    results = []

    for tool in df_clingen_ev['tool'].unique():

        # Classification concordance
        df_tool = df_clingen_master[
            df_clingen_master['tool'] == tool
        ]

        n_total = len(df_tool)

        n_concordant = len(
            df_tool[
                df_tool['classification_match'] == True
            ]
        )

        concordance = (
            n_concordant / n_total
            if n_total > 0
            else 0
        )

        # Mean Jaccard
        ev_row = df_clingen_ev[
            df_clingen_ev['tool'] == tool
        ].iloc[0]

        mean_j = ev_row['mean_criterion_jaccard']

        # Total Jaccard
        total_row = df_total_cj[
            df_total_cj['tool'] == tool
        ]

        if len(total_row) > 0:
            total_j = total_row.iloc[0][
                'raw_total_criterion_jaccard'
            ]
        else:
            total_j = None

        results.append({
            'tool': tool,
            'tool_name': TOOL_NAMES.get(tool, tool),
            'classification_concordance': concordance,
            'mean_jaccard': mean_j,
            'total_jaccard': total_j
        })

    return pd.DataFrame(results)


# =============================================================================
# PREPARE PANELS D-E
# =============================================================================

def prepare_panel_de_data(data):
    """Prepare raw vs capability-adjusted dumbbell data."""

    df_adjusted = data['adjusted']
    df_total_cj = data['total_cj']

    # Mean CJ
    df_mean = df_adjusted[
        [
            'tool',
            'raw_mean_jaccard',
            'adjusted_mean_jaccard'
        ]
    ].copy()

    df_mean['tool_name'] = df_mean['tool'].map(
        TOOL_NAMES
    )

    # Total CJ
    df_total = df_total_cj[
        [
            'tool',
            'raw_total_criterion_jaccard',
            'adjusted_total_criterion_jaccard'
        ]
    ].copy()

    df_total['tool_name'] = df_total['tool'].map(
        TOOL_NAMES
    )

    # Sort by raw Mean CJ
    df_mean = (
        df_mean
        .sort_values(
            'raw_mean_jaccard',
            ascending=False
        )
        .reset_index(drop=True)
    )

    tool_order = df_mean['tool'].tolist()

    # Apply identical order to Total CJ
    df_total['tool_order'] = df_total['tool'].map(
        {t: i for i, t in enumerate(tool_order)}
    )

    df_total = (
        df_total
        .sort_values('tool_order')
        .reset_index(drop=True)
    )

    return df_mean, df_total, tool_order


# =============================================================================
# SCATTER PANELS A-C
# =============================================================================

def plot_scatter_panel(
    ax,
    df,
    x_col,
    y_col,
    x_label,
    y_label,
    panel_label,
    show_legend=False
):

    """
    Scatter panels A-C.

    Highlighted tools:
        point size = 80

    Other tools:
        point size = 70
    """

    # -------------------------------------------------------------------------
    # Label offsets
    # -------------------------------------------------------------------------

    if panel_label == 'a':

        LABEL_OFFSETS = {
            'Genebe': (0, -10),
            'AutoGVP': (-1, 10),
            'Franklin': (0, -10),
            'TAPES': (0, -10),
            'VIP-HL': (0, 8),
            'CPSR': (-1, 10),
            'InterVar_2025': (0, -10),
        }

    elif panel_label == 'b':

        LABEL_OFFSETS = {
            'Genebe': (-8, -6),
            'AutoGVP': (-10, -7),
            'Franklin': (0, -10),
            'TAPES': (5, -9),
            'VIP-HL': (6, 3),
            'CPSR': (7, 4),
            'InterVar_2025': (4, -7),
        }

    else:

        LABEL_OFFSETS = {
            'Genebe': (-8, -6),
            'AutoGVP': (-10, -7),
            'Franklin': (-7, -7),
            'TAPES': (5, -9),
            'VIP-HL': (8, 8),
            'CPSR': (7, 4),
            'InterVar_2025': (4, -7),
        }


    # -------------------------------------------------------------------------
    # Plot points
    # -------------------------------------------------------------------------

    for _, row in df.iterrows():

        if row['tool'] in HIGHLIGHT_TOOLS:

            color = TOOL_COLORS.get(
                row['tool'],
                '#666666'
            )

            alpha = 0.95

            # Enlarged
            point_size = 68

            zorder = 4

        else:

            color = '#b8b8b8'
            alpha = 0.75

            # Enlarged
            point_size = 50

            zorder = 2


        ax.scatter(
            row[x_col],
            row[y_col],
            s=point_size,
            color=color,
            alpha=alpha,
            edgecolors='white',
            linewidth=0.5,
            zorder=zorder
        )


    # -------------------------------------------------------------------------
    # Labels for highlighted tools
    # -------------------------------------------------------------------------

    highlighted_df = df[
        df['tool'].isin(HIGHLIGHT_TOOLS)
    ]

    for _, row in highlighted_df.iterrows():

        dx, dy = LABEL_OFFSETS.get(
            row['tool'],
            (5, 5)
        )

        ax.annotate(
            row['tool_name'],
            xy=(
                row[x_col],
                row[y_col]
            ),
            xytext=(dx, dy),
            textcoords='offset points',
            fontsize=6.3,
            color='#333333',
            ha='left' if dx >= 0 else 'right',
            va='center',
            zorder=5
        )


    # -------------------------------------------------------------------------
    # Styling
    # -------------------------------------------------------------------------

    ax.set_xlabel(
        x_label,
        fontsize=8,
        fontweight='normal'
    )

    ax.set_ylabel(
        y_label,
        fontsize=8,
        fontweight='normal'
    )

    ax.tick_params(labelsize=7)

    ax.grid(
        True,
        alpha=0.2,
        linewidth=0.3
    )

    ax.set_axisbelow(True)


    # Panel label
    ax.text(
        -0.15,
        1.05,
        panel_label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight='bold',
        va='top',
        ha='left'
    )


    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)


    # Optional legend
    if show_legend:

        legend_handles = []

        for tool in sorted(HIGHLIGHT_TOOLS):

            color = TOOL_COLORS.get(
                tool,
                '#666666'
            )

            handle = plt.Line2D(
                [0],
                [0],
                marker='o',
                color='w',
                markerfacecolor=color,
                markersize=6,
                label=TOOL_NAMES[tool],
                markeredgecolor='white',
                markeredgewidth=0.5
            )

            legend_handles.append(handle)

        gray_handle = plt.Line2D(
            [0],
            [0],
            marker='o',
            color='w',
            markerfacecolor='#b8b8b8',
            markersize=6,
            label='Other tools',
            markeredgecolor='white',
            markeredgewidth=0.5
        )

        legend_handles.append(gray_handle)

        ax.legend(
            handles=legend_handles,
            loc='upper right',
            fontsize=5.5,
            frameon=True,
            edgecolor='gray',
            framealpha=0.9,
            ncol=2
        )


# =============================================================================
# DUMBBELL PANELS D-E
# =============================================================================

def plot_dumbbell_panel(
    ax,
    df,
    raw_col,
    adj_col,
    x_label,
    panel_label,
    capability_counts,
    show_legend=False
):

    """
    Raw vs capability-adjusted evidence concordance.

    Right-side annotation shows:

        supported criteria / 28

    Example:

        18/28

    Franklin is displayed as "–" because capability-adjusted
    evidence concordance was not calculated for Franklin.
    """

    n_tools = len(df)

    y_positions = np.arange(n_tools)


    # -------------------------------------------------------------------------
    # Franklin
    # -------------------------------------------------------------------------

    franklin_mask = (
        df['tool_name'] == 'Franklin'
    )

    non_franklin_mask = ~franklin_mask


    # -------------------------------------------------------------------------
    # Connector lines
    # -------------------------------------------------------------------------

    for i, (_, row) in enumerate(df.iterrows()):

        if row['tool_name'] != 'Franklin':

            raw_val = row[raw_col]
            adj_val = row[adj_col]

            ax.plot(
                [raw_val, adj_val],
                [i, i],
                color=SEMANTIC_COLORS['connector'],
                linewidth=1.5,
                zorder=1
            )


    # -------------------------------------------------------------------------
    # Raw points
    # -------------------------------------------------------------------------

    ax.scatter(
        df[raw_col],
        y_positions,
        color=SEMANTIC_COLORS['raw'],
        s=50,
        alpha=0.9,
        edgecolor='white',
        linewidth=0.5,
        zorder=3,
        label='Raw'
    )


    # -------------------------------------------------------------------------
    # Capability-adjusted points
    # -------------------------------------------------------------------------

    ax.scatter(
        df.loc[
            non_franklin_mask,
            adj_col
        ],
        y_positions[
            non_franklin_mask
        ],
        color=SEMANTIC_COLORS['adjusted'],
        s=50,
        alpha=0.9,
        edgecolor='white',
        linewidth=0.5,
        zorder=3,
        label='Capability-adjusted'
    )


    # -------------------------------------------------------------------------
    # Y-axis labels
    # -------------------------------------------------------------------------

    ax.set_yticks(y_positions)

    ax.set_yticklabels(
        df['tool_name'],
        fontsize=7
    )


    # -------------------------------------------------------------------------
    # X-axis
    # -------------------------------------------------------------------------

    ax.set_xlabel(
        x_label,
        fontsize=7,
        fontweight='normal'
    )

    ax.set_xlim(
        -0.05,
        1.05
    )

    ax.tick_params(
        axis='x',
        labelsize=7
    )

    ax.tick_params(
        axis='y',
        length=0
    )


    # -------------------------------------------------------------------------
    # Grid
    # -------------------------------------------------------------------------

    ax.grid(
        True,
        axis='x',
        alpha=0.2,
        linewidth=0.3
    )

    ax.set_axisbelow(True)

    ax.spines['left'].set_visible(True)


    # -------------------------------------------------------------------------
    # Panel label
    # -------------------------------------------------------------------------

    ax.text(
        -0.12,
        1.05,
        panel_label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight='bold',
        va='top',
        ha='left'
    )


    # =========================================================================
    # SUPPORTED CRITERIA COLUMN
    # =========================================================================

    # Position outside the plotting area using axes coordinates.
    # 1.00 = right edge of plotting area.
    # 1.12 = 12% beyond the right edge.
    criteria_x = 1.12

    for i, (_, row) in enumerate(df.iterrows()):

        tool = row['tool']

        # Franklin was intentionally excluded from capability adjustment
        if tool == 'Franklin':
            capability_text = '–'
        else:
            n_supported = capability_counts.get(tool)

            if n_supported is None:
                capability_text = '–'
            else:
                capability_text = f'{int(n_supported)}/28'

        # Convert data-coordinate row position to axes-coordinate y position
        y_axes = ax.transAxes.inverted().transform(
            ax.transData.transform((0, i))
        )[1]

        ax.text(
            criteria_x,
            y_axes,
            capability_text,
            transform=ax.transAxes,
            ha='center',
            va='center',
            fontsize=6.5,
            color='#444444',
            clip_on=False
        )

    # Column heading
    ax.text(
        criteria_x,
        1.02,
        'Supported \ncriteria',
        transform=ax.transAxes,
        ha='center',
        va='bottom',
        fontsize=6.5,
        fontweight='bold',
        color='#333333',
        linespacing=1.0,
        clip_on=False
    )

    # -------------------------------------------------------------------------
    # Legend
    # -------------------------------------------------------------------------

    if show_legend:
        ax.legend(
            loc='upper right',
            fontsize=6,
            frameon=True,
            edgecolor='gray',
            framealpha=0.9
        )


# =============================================================================
# GENERATE FIGURE
# =============================================================================

def generate_five_panel_figure_v3(
    data,
    output_dir,
    capability_counts
):

    print("\n" + "=" * 80)
    print("GENERATING FIVE-PANEL FIGURE 5 - V3")
    print("=" * 80)


    # -------------------------------------------------------------------------
    # Prepare data
    # -------------------------------------------------------------------------

    print("\nPreparing data...")

    df_scatter = prepare_panel_abc_data(
        data
    )

    df_mean, df_total, tool_order = (
        prepare_panel_de_data(data)
    )

    print(
        f"  ✓ Panel A-C data: "
        f"{len(df_scatter)} tools"
    )

    print(
        f"  ✓ Panel D data: "
        f"{len(df_mean)} tools"
    )

    print(
        f"  ✓ Panel E data: "
        f"{len(df_total)} tools"
    )


    # -------------------------------------------------------------------------
    # Figure
    # -------------------------------------------------------------------------

    fig = plt.figure(
        figsize=(13, 8.5)
    )


    outer = GridSpec(
        2,
        1,
        figure=fig,
        height_ratios=[1, 1],
        hspace=0.42,
        left=0.08,
        right=0.96,
        top=0.96,
        bottom=0.08
    )


    # -------------------------------------------------------------------------
    # Top row
    # -------------------------------------------------------------------------

    top = outer[0].subgridspec(
        1,
        3,
        wspace=0.18
    )


    ax_a = fig.add_subplot(
        top[0, 0]
    )

    ax_b = fig.add_subplot(
        top[0, 1]
    )

    ax_c = fig.add_subplot(
        top[0, 2]
    )


    # -------------------------------------------------------------------------
    # Bottom row
    # -------------------------------------------------------------------------

    bottom = outer[1].subgridspec(
        1,
        5,
        width_ratios=[0.45, 2.0, 1.15, 2.0, 0.75],
        wspace=0.05
    )


    ax_d = fig.add_subplot(
        bottom[0, 1]
    )

    ax_e = fig.add_subplot(
        bottom[0, 3]
    )


    print("\nPlotting panels...")


    # =========================================================================
    # PANEL A
    # =========================================================================

    print(
        "  ✓ Panel A: "
        "Mean CJ vs Total CJ"
    )

    plot_scatter_panel(
        ax_a,
        df_scatter,
        'mean_jaccard',
        'total_jaccard',
        'Mean criterion Jaccard',
        'Total criterion Jaccard',
        'a'
    )


    # =========================================================================
    # PANEL B
    # =========================================================================

    print(
        "  ✓ Panel B: "
        "Classification concordance vs Mean CJ"
    )

    plot_scatter_panel(
        ax_b,
        df_scatter,
        'classification_concordance',
        'mean_jaccard',
        'Classification concordance (valid predictions)',
        'Mean criterion Jaccard',
        'b'
    )


    # =========================================================================
    # PANEL C
    # =========================================================================

    print(
        "  ✓ Panel C: "
        "Classification concordance vs Total CJ"
    )

    plot_scatter_panel(
        ax_c,
        df_scatter,
        'classification_concordance',
        'total_jaccard',
        'Classification concordance (valid predictions)',
        'Total criterion Jaccard',
        'c',
        show_legend=False
    )


    # =========================================================================
    # PANEL D
    # =========================================================================

    print(
        "  ✓ Panel D: "
        "Raw vs Adjusted Mean CJ"
    )

    plot_dumbbell_panel(
        ax_d,
        df_mean,
        'raw_mean_jaccard',
        'adjusted_mean_jaccard',
        'Mean criterion Jaccard',
        'd',
        capability_counts,
        show_legend=True
    )


    # =========================================================================
    # PANEL E
    # =========================================================================

    print(
        "  ✓ Panel E: "
        "Raw vs Adjusted Total CJ"
    )

    plot_dumbbell_panel(
        ax_e,
        df_total,
        'raw_total_criterion_jaccard',
        'adjusted_total_criterion_jaccard',
        'Total criterion Jaccard',
        'e',
        capability_counts,
        show_legend=True
    )


    # -------------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------------

    print("\nSaving figure...")

    output_dir = Path(
        output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    png_path = (
        output_dir /
        'Figure5_five_panel_preview_v3.png'
    )

    pdf_path = (
        output_dir /
        'Figure5_five_panel_preview_v3.pdf'
    )


    plt.savefig(
        png_path,
        dpi=300,
        bbox_inches='tight'
    )

    plt.savefig(
        pdf_path,
        bbox_inches='tight'
    )


    print(
        f"  ✓ PNG: {png_path}"
    )

    print(
        f"  ✓ PDF: {pdf_path}"
    )


    plt.close()


    return (
        df_scatter,
        df_mean,
        df_total,
        tool_order
    )


# =============================================================================
# VALIDATION REPORT
# =============================================================================

def generate_report(
    df_scatter,
    df_mean,
    df_total,
    tool_order,
    capability_counts
):

    print("\n" + "=" * 80)
    print("FIGURE 5 V3 VALIDATION REPORT")
    print("=" * 80)


    print("\n1. OUTPUT:")
    print(
        "   ✓ Figure5_five_panel_preview_v3.png"
    )
    print(
        "   ✓ Figure5_five_panel_preview_v3.pdf"
    )


    print("\n2. PANELS A-C:")
    print(
        "   ✓ Highlighted point size = 68"
    )
    print(
        "   ✓ Context point size = 50"
    )
    print(
        "   ✓ No metric recalculation"
    )


    print("\n3. PANELS D-E:")
    print(
        "   ✓ Raw and adjusted values unchanged"
    )
    print(
        "   ✓ Supported criteria shown as n/28"
    )
    print(
        "   ✓ Franklin shown as –"
    )


    print("\n4. CAPABILITY COUNTS USED:")

    for tool in tool_order:

        display_name = TOOL_NAMES.get(
            tool,
            tool
        )

        if tool == 'Franklin':

            value = '–'

        else:

            count = capability_counts.get(
                tool,
                None
            )

            if count is None:
                value = '–'
            else:
                value = (
                    f'{count}/'
                    f'{TOTAL_ACMG_CRITERIA}'
                )

        print(
            f"   {display_name:<20} {value}"
        )


    print("\n5. DATA INTEGRITY:")
    print(
        "   ✓ No performance metrics recalculated"
    )
    print(
        "   ✓ No Jaccard values recalculated"
    )
    print(
        "   ✓ No population filters changed"
    )
    print(
        "   ✓ Existing capability-adjusted values unchanged"
    )
    print(
        "   ✓ Capability matrix used only for displayed n/28 counts"
    )


    print("\n" + "=" * 80)
    print("✓ FIGURE 5 V3 COMPLETE")
    print("=" * 80)


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 80)
    print(
        "FIGURE 5 - FIVE-PANEL DESIGN V3"
    )
    print("=" * 80)


    # -------------------------------------------------------------------------
    # Load existing figure data
    # -------------------------------------------------------------------------

    print("\nLoading figure data...")

    data = load_data()

    print(
        "  ✓ Existing figure data loaded"
    )


    # -------------------------------------------------------------------------
    # Load capability matrix
    # -------------------------------------------------------------------------

    capability_counts = (
        load_capability_counts()
    )


    # -------------------------------------------------------------------------
    # Generate figure
    # -------------------------------------------------------------------------

    (
        df_scatter,
        df_mean,
        df_total,
        tool_order
    ) = generate_five_panel_figure_v3(
        data,
        'outputs/figures',
        capability_counts
    )


    # -------------------------------------------------------------------------
    # Validation report
    # -------------------------------------------------------------------------

    generate_report(
        df_scatter,
        df_mean,
        df_total,
        tool_order,
        capability_counts
    )


if __name__ == '__main__':
    main()
