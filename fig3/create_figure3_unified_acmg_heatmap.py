#!/usr/bin/env python3
"""
Unified ACMG Criteria Heatmap - Two Panel Design
Panel A: ClinGen | Panel B: FOXL2
All criteria in single Y-axis with category separators
No dendrograms, clean design matching benchmark figure
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path

# Global font settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['text.color'] = '#222222'
plt.rcParams['axes.labelcolor'] = '#222222'
plt.rcParams['font.size'] = 8

print("=" * 80)
print("UNIFIED ACMG CRITERIA HEATMAP - TWO PANEL DESIGN")
print("=" * 80)

# Base directory and output - adjusted for fig2 folder structure
script_dir = Path(__file__).parent
base_dir = script_dir.parent.parent / 'analysis_last'
# Output to analysis_scripts/figures directory
output_dir = script_dir.parent / 'figures'
output_dir.mkdir(parents=True, exist_ok=True)

# ACMG Categories (only the 3 main ones for unified display)
CATEGORIES = {
    'Population': ['BA1', 'BS1', 'BS2', 'PM2', 'PS4'],
    'Computational and Predictive': ['BP1', 'BP3', 'BP4', 'BP7', 'PP3', 'PM4', 'PM5', 'PS1', 'PVS1'],
    'Functional': ['BS3', 'PP2', 'PM1', 'PS3']
}

# All tools
ALL_TOOLS = ['Genebe', 'BIAS', 'Franklin', 'CharGer_Local', 'CharGer_Online',
             'InterVar_2025', 'InterVar_2018', 'AutoGVP', 'DiabloACMG', 
             'TAPES', 'Exomiser', 'CPSR', 'VIP-HL', 'CancerSIGVAR']


def load_jaccard_data(dataset_name):
    """Load Jaccard similarity data"""
    # First try local directory (for standalone use)
    if dataset_name == 'ClinGen':
        local_file = script_dir / 'jaccard_per_criterion_clingen_28012026.csv'
        original_file = base_dir / 'results' / 'jaccard_per_criterion_clingen_28012026.csv'
    elif dataset_name == 'FOXL2':
        local_file = script_dir / 'jaccard_per_criterion_foxl2.csv'
        original_file = base_dir / 'results' / 'jaccard_per_criterion_foxl2.csv'
    else:
        raise ValueError(f"Unknown dataset: {dataset_name}")
    
    # Use local file if exists, otherwise use original location
    if local_file.exists():
        file_path = local_file
    elif original_file.exists():
        file_path = original_file
    else:
        print(f"   WARNING: Data file not found in {local_file} or {original_file}")
        return None
    
    return pd.read_csv(file_path)


def create_jaccard_matrix(df, criteria, tools):
    """Create Jaccard similarity matrix"""
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


def create_unified_heatmap():
    """Create unified two-panel heatmap"""
    print("\n1. Loading data...")
    
    # Load data for both datasets
    clingen_df = load_jaccard_data('ClinGen')
    foxl2_df = load_jaccard_data('FOXL2')
    
    if clingen_df is None or foxl2_df is None:
        print("   ERROR: Could not load data")
        return
    
    print("\n2. Preparing matrices...")
    
    # Combine all criteria in order
    all_criteria = []
    category_boundaries = []
    category_labels = []
    current_pos = 0
    
    for cat_name, criteria_list in CATEGORIES.items():
        all_criteria.extend(criteria_list)
        category_boundaries.append(current_pos)
        category_labels.append(cat_name)
        current_pos += len(criteria_list)
    
    category_boundaries.append(current_pos)  # Final boundary
    
    print(f"   Total criteria: {len(all_criteria)}")
    print(f"   Categories: {list(CATEGORIES.keys())}")
    
    # Create matrices for both datasets
    clingen_matrix, clingen_tools = create_jaccard_matrix(clingen_df, all_criteria, ALL_TOOLS)
    foxl2_matrix, foxl2_tools = create_jaccard_matrix(foxl2_df, all_criteria, ALL_TOOLS)
    
    print(f"   ClinGen: {clingen_matrix.shape}")
    print(f"   FOXL2: {foxl2_matrix.shape}")
    
    print("\n3. Creating figure...")
    
    # Create figure with two panels side by side
    fig = plt.figure(figsize=(16, 10))
    
    # Create grid: 2 columns for the two panels
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 1], wspace=0.3)
    
    # Panel A: ClinGen
    ax_a = fig.add_subplot(gs[0])
    
    # Panel B: FOXL2
    ax_b = fig.add_subplot(gs[1])
    
    # Color palette matching benchmark figure (RdBu_r)
    cmap = 'RdBu_r'
    
    # Plot Panel A (ClinGen)
    print("   Plotting Panel A (ClinGen)...")
    sns.heatmap(
        clingen_matrix,
        ax=ax_a,
        cmap=cmap,
        vmin=0,
        vmax=1,
        cbar=False,  # We'll add a single colorbar later
        xticklabels=clingen_tools,
        yticklabels=all_criteria,
        linewidths=0.5,
        linecolor='white',
        square=False,
        annot=True,
        fmt='.2f',
        annot_kws={'fontsize': 6}
    )
    
    ax_a.set_title('ClinGen Master Cohort (n=11,409)', 
                   fontsize=12, fontweight='bold', pad=10)
    ax_a.set_xlabel('Tools', fontsize=10, fontweight='bold')
    ax_a.set_ylabel('')  # Remove ACMG Criterion label
    ax_a.set_xticklabels(ax_a.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax_a.set_yticklabels(ax_a.get_yticklabels(), rotation=0, fontsize=8)
    
    # Add panel label 'a' at top-left (like Figure1)
    ax_a.text(-0.15, 1.05, 'a', transform=ax_a.transAxes,
              fontsize=20, fontweight='bold', va='top', ha='right')
    
    # Plot Panel B (FOXL2)
    print("   Plotting Panel B (FOXL2)...")
    im = sns.heatmap(
        foxl2_matrix,
        ax=ax_b,
        cmap=cmap,
        vmin=0,
        vmax=1,
        cbar=True,
        cbar_kws={'label': 'Jaccard Similarity Coefficient', 'shrink': 0.8},
        xticklabels=foxl2_tools,
        yticklabels=all_criteria,
        linewidths=0.5,
        linecolor='white',
        square=False,
        annot=True,
        fmt='.2f',
        annot_kws={'fontsize': 6}
    )
    
    ax_b.set_title('FOXL2 Validation Cohort (n=288)', 
                   fontsize=12, fontweight='bold', pad=10)
    ax_b.set_xlabel('Tools', fontsize=10, fontweight='bold')
    ax_b.set_ylabel('')  # No Y label for second panel
    ax_b.set_xticklabels(ax_b.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax_b.set_yticklabels(ax_b.get_yticklabels(), rotation=0, fontsize=8)
    
    # Add panel label 'b' at top-left (like Figure1)
    ax_b.text(-0.15, 1.05, 'b', transform=ax_b.transAxes,
              fontsize=20, fontweight='bold', va='top', ha='right')
    
    # Add category separators and labels to both panels
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
    
    # No overall title (removed as requested)
    
    # Save
    print("\n4. Saving figure...")
    output_file = output_dir / 'Figure3_Unified_ACMG_Heatmap_TwoPanel.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"   ✓ Saved: {output_file}")
    
    output_file_pdf = output_dir / 'Figure3_Unified_ACMG_Heatmap_TwoPanel.pdf'
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight')
    print(f"   ✓ Saved: {output_file_pdf}")
    
    plt.close()


# Main execution
if __name__ == '__main__':
    create_unified_heatmap()
    
    print("\n" + "=" * 80)
    print("✓ UNIFIED HEATMAP CREATED")
    print(f"✓ Output directory: {output_dir}")
    print("=" * 80)
