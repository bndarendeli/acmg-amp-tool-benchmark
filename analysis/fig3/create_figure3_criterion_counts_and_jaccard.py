#!/usr/bin/env python3
"""
Figure 3: ACMG Criterion Counts and Jaccard Similarity

Two-panel design for both ClinGen and FOXL2:
- Panel a: Criterion assignment counts (tools + Ground truth column)
- Panel b: Criterion-level Jaccard similarity

Uses corrected pipeline with proper missing-reference handling.
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path
import sys

# Import criterion parser
sys.path.insert(0, str(Path(__file__).parent.parent / 'fig2'))
from utils_criteria_parser import parse_criteria_string

# Global settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 7

print("="*80)
print("FIGURE 3: CRITERION COUNTS AND JACCARD SIMILARITY")
print("="*80)

# Paths
script_dir = Path(__file__).parent
MASTER_TABLE = script_dir.parent / 'fig5' / 'outputs' / 'fig4_master_variant_tool.csv'
JACCARD_CLINGEN = script_dir / 'jaccard_per_criterion_clingen_28012026.csv'
JACCARD_FOXL2 = script_dir / 'jaccard_per_criterion_foxl2.csv'
output_dir = script_dir / 'figures'
output_dir.mkdir(parents=True, exist_ok=True)

# Selected ACMG/AMP criteria (from original Figure 3 design)
# Organized by category
CATEGORIES = {
    'Population': ['BA1', 'BS1', 'BS2', 'PM2', 'PS4'],
    'Computational and Predictive': ['BP1', 'BP3', 'BP4', 'BP7', 'PP3', 'PM4', 'PM5', 'PS1', 'PVS1'],
    'Functional': ['BS3', 'PP2', 'PM1', 'PS3']
}

# Flatten to get ordered list
SELECTED_CRITERIA = []
for cat_name, criteria_list in CATEGORIES.items():
    SELECTED_CRITERIA.extend(criteria_list)

# All tools
ALL_TOOLS = [
    'InterVar_2018', 'InterVar_2025', 'BIAS', 'CharGer_Local', 'CharGer_Online',
    'DiabloACMG', 'Exomiser', 'Genebe', 'TAPES', 'Franklin', 'AutoGVP', 'VIP-HL',
    'CancerSIGVAR', 'CPSR'
]


def is_criteria_available(criteria_str):
    """Check if criteria string represents available data (not missing)"""
    if pd.isna(criteria_str):
        return False
    criteria_str = str(criteria_str).strip()
    if criteria_str in ['', '.', 'Unknown', 'Not_Provided', 'None', 'nan', 'NA']:
        return False
    return True


def calculate_criterion_counts(df, dataset_name, tools):
    """
    Calculate how many times each criterion was assigned by each tool.
    
    Returns:
        count_matrix: tools × criteria matrix of counts
        ground_truth_counts: criteria counts from reference
        available_tools: tools that have data for this dataset
    """
    print(f"\nCalculating criterion counts for {dataset_name}...")
    
    # Filter to dataset
    df_dataset = df[df['dataset'] == dataset_name].copy()
    print(f"  Total variant-tool pairs: {len(df_dataset)}")
    
    # Initialize count matrices
    count_matrix = np.zeros((len(tools), len(SELECTED_CRITERIA)), dtype=int)
    ground_truth_counts = np.zeros(len(SELECTED_CRITERIA), dtype=int)
    
    # Track which tools have data
    available_tools = []
    tool_indices = {}
    
    for tool_idx, tool in enumerate(tools):
        df_tool = df_dataset[df_dataset['tool'] == tool]
        if len(df_tool) == 0:
            continue
        
        available_tools.append(tool)
        tool_indices[tool] = len(available_tools) - 1
    
    # Count criterion assignments
    for _, row in df_dataset.iterrows():
        tool = row['tool']
        if tool not in tool_indices:
            continue
        
        tool_idx = tool_indices[tool]
        
        # Parse tool criteria
        tool_criteria_str = row['tool_criteria']
        if is_criteria_available(tool_criteria_str):
            tool_criteria = parse_criteria_string(tool_criteria_str)
            for criterion in tool_criteria:
                if criterion in SELECTED_CRITERIA:
                    crit_idx = SELECTED_CRITERIA.index(criterion)
                    count_matrix[tool_idx, crit_idx] += 1
    
    # Count ground truth (reference) criteria
    # Use unique variants only
    df_variants = df_dataset[['variant_id', 'reference_criteria']].drop_duplicates()
    
    for _, row in df_variants.iterrows():
        ref_criteria_str = row['reference_criteria']
        if is_criteria_available(ref_criteria_str):
            ref_criteria = parse_criteria_string(ref_criteria_str)
            for criterion in ref_criteria:
                if criterion in SELECTED_CRITERIA:
                    crit_idx = SELECTED_CRITERIA.index(criterion)
                    ground_truth_counts[crit_idx] += 1
    
    # Trim count matrix to available tools
    count_matrix = count_matrix[:len(available_tools), :]
    
    print(f"  Available tools: {len(available_tools)}")
    print(f"  Count matrix shape: {count_matrix.shape}")
    print(f"  Ground truth counts: {ground_truth_counts.sum()} total assignments")
    
    return count_matrix, ground_truth_counts, available_tools


def load_jaccard_matrix(filepath, criteria, tools):
    """Load Jaccard similarity matrix from CSV"""
    if not filepath.exists():
        print(f"  WARNING: Jaccard file not found: {filepath}")
        return None, []
    
    df = pd.read_csv(filepath)
    
    matrix = np.zeros((len(criteria), len(tools)))
    available_tools = []
    
    for j, tool in enumerate(tools):
        tool_data = df[df['Tool'] == tool]
        if len(tool_data) > 0:
            available_tools.append(tool)
            for i, criterion in enumerate(criteria):
                value = tool_data[tool_data['Criterion'] == criterion]['Jaccard_Similarity']
                if len(value) > 0:
                    matrix[i, len(available_tools)-1] = value.values[0]
    
    # Trim matrix to available tools
    matrix = matrix[:, :len(available_tools)]
    
    return matrix, available_tools


def create_figure_for_dataset(dataset_name, count_matrix, ground_truth_counts, 
                               jaccard_matrix, tools, criteria, output_prefix):
    """
    Create two-panel figure for a dataset:
    - Panel a: Criterion counts (tools + Ground truth)
    - Panel b: Jaccard similarity
    """
    print(f"\nCreating figure for {dataset_name}...")
    
    # Calculate category boundaries
    category_boundaries = []
    category_labels = []
    current_pos = 0
    
    for cat_name, criteria_list in CATEGORIES.items():
        category_boundaries.append(current_pos)
        category_labels.append(cat_name)
        current_pos += len(criteria_list)
    category_boundaries.append(current_pos)  # Final boundary
    
    # Create figure matching old layout
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1.1, 1], wspace=0.3)
    
    # Panel a: Criterion counts
    ax_a = fig.add_subplot(gs[0])
    
    # Combine count matrix with ground truth column
    # Transpose so criteria are on Y-axis
    count_matrix_t = count_matrix.T  # criteria × tools
    
    # Add ground truth as last column
    count_with_gt = np.column_stack([count_matrix_t, ground_truth_counts])
    
    # Create labels
    x_labels_a = tools + ['Ground truth']
    
    # Plot count heatmap (Blues color scale)
    sns.heatmap(
        count_with_gt,
        ax=ax_a,
        cmap='Blues',  # Light to dark blue
        cbar=True,
        cbar_kws={'label': 'Criterion assignment count', 'shrink': 0.8, 'pad': 0.02},
        xticklabels=x_labels_a,
        yticklabels=criteria,
        linewidths=0.5,
        linecolor='white',
        square=False,
        annot=True,
        fmt='d',
        annot_kws={'fontsize': 6}
    )
    
    # Add vertical separator before Ground truth column (thinner)
    n_tools = len(tools)
    ax_a.axvline(x=n_tools, color='black', linewidth=1.5, linestyle='-')
    
    # Adjust annotation text colors for readability (black on light, white on dark)
    # Also hide text for cells with count = 0
    count_max = count_with_gt.max()
    count_min = count_with_gt.min()
    threshold = (count_max - count_min) * 0.5 + count_min
    
    for text in ax_a.texts:
        # Get the text value and position
        try:
            value = float(text.get_text())
            # Hide text for zero counts
            if value == 0:
                text.set_text('')
            elif value > threshold:
                text.set_color('white')
            else:
                text.set_color('black')
        except:
            pass
    
    ax_a.set_xlabel('', fontsize=8)
    ax_a.set_ylabel('', fontsize=8)
    ax_a.set_xticklabels(ax_a.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax_a.set_yticklabels(ax_a.get_yticklabels(), rotation=0, fontsize=8)
    
    # Panel label (matching old style)
    ax_a.text(-0.15, 1.05, 'a', transform=ax_a.transAxes,
              fontsize=20, fontweight='bold', va='top', ha='right')
    
    # Panel b: Jaccard similarity
    ax_b = fig.add_subplot(gs[1])
    
    # Plot Jaccard heatmap (original RdBu_r color scale - keep as original)
    sns.heatmap(
        jaccard_matrix,
        ax=ax_b,
        cmap='RdBu_r',  # Original red-blue color scale
        vmin=0,
        vmax=1,
        cbar=True,
        cbar_kws={'label': 'Criterion-level Jaccard similarity', 'shrink': 0.8},
        xticklabels=tools,
        yticklabels=criteria,
        linewidths=0.5,
        linecolor='white',
        square=False,
        annot=True,
        fmt='.2f',
        annot_kws={'fontsize': 6}
    )
    
    # Adjust annotation text colors for readability in Jaccard panel
    # For RdBu_r: values near 0.5 are lightest (white), extremes are dark
    for text in ax_b.texts:
        try:
            value = float(text.get_text())
            # Light colors (near 0.5) get black text, dark colors get white text
            if 0.3 < value < 0.7:
                text.set_color('black')
            else:
                text.set_color('white')
        except:
            pass
    
    ax_b.set_xlabel('', fontsize=8)
    ax_b.set_ylabel('', fontsize=8)
    ax_b.set_xticklabels(ax_b.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax_b.set_yticklabels(ax_b.get_yticklabels(), rotation=0, fontsize=8)
    
    # Panel label (matching old style)
    ax_b.text(-0.15, 1.05, 'b', transform=ax_b.transAxes,
              fontsize=20, fontweight='bold', va='top', ha='right')
    
    # Add category separators and labels to both panels (matching old figure)
    for ax in [ax_a, ax_b]:
        # Add horizontal lines to separate categories
        for boundary in category_boundaries[1:-1]:  # Skip first and last
            ax.axhline(y=boundary, color='black', linewidth=2, linestyle='-')
        
        # Add thin brackets to highlight category regions
        for i, (cat_name, start, end) in enumerate(zip(category_labels, 
                                                        category_boundaries[:-1], 
                                                        category_boundaries[1:])):
            # Draw thin bracket on the left side
            bracket_x = -1.8
            bracket_width = 0.15
            
            # Vertical line (bracket spine)
            ax.plot([bracket_x, bracket_x], [start, end], 
                   color='black', linewidth=1.5, clip_on=False)
            
            # Top horizontal line (bracket cap)
            ax.plot([bracket_x, bracket_x + bracket_width], [start, start], 
                   color='black', linewidth=1.5, clip_on=False)
            
            # Bottom horizontal line (bracket cap)
            ax.plot([bracket_x, bracket_x + bracket_width], [end, end], 
                   color='black', linewidth=1.5, clip_on=False)
        
        # Add category labels on the left side (further from y-axis labels)
        for i, (cat_name, start, end) in enumerate(zip(category_labels, 
                                                        category_boundaries[:-1], 
                                                        category_boundaries[1:])):
            mid_point = (start + end) / 2
            # Add text further to the left to avoid overlap with criterion labels
            ax.text(-2.8, mid_point, cat_name, 
                   rotation=90, va='center', ha='right',
                   fontsize=9, fontweight='bold')
    
    # Save
    png_path = output_dir / f"{output_prefix}.png"
    pdf_path = output_dir / f"{output_prefix}.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {png_path.name}")
    
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {pdf_path.name}")
    
    plt.close()
    
    return count_with_gt, jaccard_matrix


def save_matrices_as_csv(dataset_name, count_matrix, ground_truth_counts, 
                         jaccard_matrix, tools, criteria):
    """Save count and Jaccard matrices as CSV files"""
    print(f"\nSaving matrices for {dataset_name}...")
    
    # Count matrix with ground truth
    count_with_gt = np.column_stack([count_matrix.T, ground_truth_counts])
    df_counts = pd.DataFrame(
        count_with_gt,
        index=criteria,
        columns=tools + ['Ground_truth']
    )
    df_counts.index.name = 'Criterion'
    
    count_path = output_dir / f"figure3_{dataset_name.lower()}_criterion_counts.csv"
    df_counts.to_csv(count_path)
    print(f"  ✓ Saved: {count_path.name}")
    
    # Jaccard matrix
    df_jaccard = pd.DataFrame(
        jaccard_matrix,
        index=criteria,
        columns=tools
    )
    df_jaccard.index.name = 'Criterion'
    
    jaccard_path = output_dir / f"figure3_{dataset_name.lower()}_jaccard.csv"
    df_jaccard.to_csv(jaccard_path)
    print(f"  ✓ Saved: {jaccard_path.name}")


def main():
    """Main execution"""
    
    # Load master table
    print("\nLoading master table...")
    if not MASTER_TABLE.exists():
        print(f"  ERROR: Master table not found: {MASTER_TABLE}")
        return
    
    df = pd.read_csv(MASTER_TABLE)
    print(f"  Loaded {len(df)} variant-tool pairs")
    
    # ========================================================================
    # ClinGen Main Figure
    # ========================================================================
    print("\n" + "="*80)
    print("CLINGEN MAIN FIGURE 3")
    print("="*80)
    
    # Calculate counts
    clingen_counts, clingen_gt, clingen_tools = calculate_criterion_counts(
        df, 'clingen_28012026', ALL_TOOLS
    )
    
    # Load Jaccard
    clingen_jaccard, clingen_jaccard_tools = load_jaccard_matrix(
        JACCARD_CLINGEN, SELECTED_CRITERIA, ALL_TOOLS
    )
    
    if clingen_jaccard is None:
        print("  ERROR: Could not load ClinGen Jaccard data")
        return
    
    # Verify tool consistency
    if clingen_tools != clingen_jaccard_tools:
        print(f"  WARNING: Tool mismatch between counts and Jaccard")
        print(f"    Counts: {clingen_tools}")
        print(f"    Jaccard: {clingen_jaccard_tools}")
        # Use intersection
        common_tools = [t for t in clingen_tools if t in clingen_jaccard_tools]
        print(f"    Using common tools: {len(common_tools)}")
        
        # Reorder matrices
        count_indices = [clingen_tools.index(t) for t in common_tools]
        jaccard_indices = [clingen_jaccard_tools.index(t) for t in common_tools]
        
        clingen_counts = clingen_counts[count_indices, :]
        clingen_jaccard = clingen_jaccard[:, jaccard_indices]
        clingen_tools = common_tools
    
    # Create figure
    create_figure_for_dataset(
        'ClinGen',
        clingen_counts,
        clingen_gt,
        clingen_jaccard,
        clingen_tools,
        SELECTED_CRITERIA,
        'Figure3_ClinGen_Counts_Jaccard'
    )
    
    # Save matrices
    save_matrices_as_csv(
        'ClinGen',
        clingen_counts,
        clingen_gt,
        clingen_jaccard,
        clingen_tools,
        SELECTED_CRITERIA
    )
    
    # ========================================================================
    # FOXL2 Supplementary Figure
    # ========================================================================
    print("\n" + "="*80)
    print("FOXL2 SUPPLEMENTARY FIGURE 3")
    print("="*80)
    
    # Calculate counts
    foxl2_counts, foxl2_gt, foxl2_tools = calculate_criterion_counts(
        df, 'foxl2', ALL_TOOLS
    )
    
    # Load Jaccard
    foxl2_jaccard, foxl2_jaccard_tools = load_jaccard_matrix(
        JACCARD_FOXL2, SELECTED_CRITERIA, ALL_TOOLS
    )
    
    if foxl2_jaccard is None:
        print("  ERROR: Could not load FOXL2 Jaccard data")
        return
    
    # Verify tool consistency
    if foxl2_tools != foxl2_jaccard_tools:
        print(f"  WARNING: Tool mismatch between counts and Jaccard")
        print(f"    Counts: {foxl2_tools}")
        print(f"    Jaccard: {foxl2_jaccard_tools}")
        # Use intersection
        common_tools = [t for t in foxl2_tools if t in foxl2_jaccard_tools]
        print(f"    Using common tools: {len(common_tools)}")
        
        # Reorder matrices
        count_indices = [foxl2_tools.index(t) for t in common_tools]
        jaccard_indices = [foxl2_jaccard_tools.index(t) for t in common_tools]
        
        foxl2_counts = foxl2_counts[count_indices, :]
        foxl2_jaccard = foxl2_jaccard[:, jaccard_indices]
        foxl2_tools = common_tools
    
    # Create figure
    create_figure_for_dataset(
        'FOXL2',
        foxl2_counts,
        foxl2_gt,
        foxl2_jaccard,
        foxl2_tools,
        SELECTED_CRITERIA,
        'Figure3_FOXL2_Counts_Jaccard'
    )
    
    # Save matrices
    save_matrices_as_csv(
        'FOXL2',
        foxl2_counts,
        foxl2_gt,
        foxl2_jaccard,
        foxl2_tools,
        SELECTED_CRITERIA
    )
    
    print("\n" + "="*80)
    print("✓ FIGURE 3 GENERATION COMPLETE")
    print("="*80)
    print(f"\nOutput directory: {output_dir}")
    print("\nGenerated files:")
    print("  ClinGen:")
    print("    - Figure3_ClinGen_Counts_Jaccard.png")
    print("    - Figure3_ClinGen_Counts_Jaccard.pdf")
    print("    - figure3_clingen_criterion_counts.csv")
    print("    - figure3_clingen_jaccard.csv")
    print("  FOXL2:")
    print("    - Figure3_FOXL2_Counts_Jaccard.png")
    print("    - Figure3_FOXL2_Counts_Jaccard.pdf")
    print("    - figure3_foxl2_criterion_counts.csv")
    print("    - figure3_foxl2_jaccard.csv")
    print("\n" + "="*80)


if __name__ == '__main__':
    main()
