#!/usr/bin/env python3
"""
Figure 4: Classification Performance - DOT PLOT DESIGN
2x2 Multi-Panel Layout with Horizontal Multi-Marker Dot Plots

Design:
- Weighted F1 (3-class): Large filled circle (primary metric)
- P/LP Recall: Smaller marker
- B/LB Recall: Smaller marker  
- VUS Recall: Smaller marker
- Call Rate removed (covered in Figure 2)
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from pathlib import Path

# Typography
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 9
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 11
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9

print("=" * 80)
print("FIGURE 4: CLASSIFICATION PERFORMANCE (DOT PLOT DESIGN)")
print("=" * 80)

# Paths
base_dir = Path(__file__).parent
data_dir = base_dir / 'data' / 'corrected'
output_dir = base_dir / 'figures'
output_dir.mkdir(parents=True, exist_ok=True)

# Load corrected metrics
print("\n1. Loading corrected metrics...")

cohorts = {
    'ClinGen': {
        'n': 11409,
        'file': 'ClinGen_corrected_metrics.csv',
        'panel_label': 'a',
        'title': 'ClinGen Master Cohort'
    },
    'FOXL2': {
        'n': 288,
        'file': 'FOXL2_corrected_metrics.csv',
        'panel_label': 'b',
        'title': 'FOXL2 Disease-Specific Cohort'
    },
    'HGMD–ClinVar HL': {
        'n': 1948,
        'file': 'HGMD+ClinVar_HL_corrected_metrics.csv',
        'panel_label': 'c',
        'title': 'HGMD–ClinVar Hearing Loss'
    },
    'HGMD–ClinVar Cancer': {
        'n': 1788,
        'file': 'HGMD+ClinVar_Cancer_corrected_metrics.csv',
        'panel_label': 'd',
        'title': 'HGMD–ClinVar Cancer Predisposition'
    }
}

cohort_data = {}
for cohort_name, cohort_info in cohorts.items():
    file_path = data_dir / cohort_info['file']
    if file_path.exists():
        df = pd.read_csv(file_path)
        
        # Remove artifact rows
        df = df[~df['Tool'].isin(['analysis', 'Genebe_Gene'])]
        
        # Sort by Weighted F1 (descending for top-to-bottom display)
        df = df.sort_values('Weighted F1 (3-class)', ascending=False)
        
        cohort_data[cohort_name] = df
        print(f"   ✓ {cohort_name}: {len(df)} tools, mean F1={df['Weighted F1 (3-class)'].mean():.3f}")
    else:
        print(f"   ✗ {cohort_name}: File not found - {file_path}")

if not cohort_data:
    print("\n✗ No data files found. Please run recalculate_correct_metrics.py first.")
    exit(1)

# Refined publication-quality palette with clear visual hierarchy
metric_styles = {
    'Weighted F1 (3-class)': {
        'color': '#4B0082',      # Deep indigo/purple - primary metric
        'marker': 'o',           # Circle
        'size': 75,              # Clearly larger but not oversized
        'edgewidth': 0.5,
        'label': 'Weighted F1 (3-class)',
        'zorder': 4
    },
    'P/LP Recall': {
        'color': '#D94F30',      # Muted vermillion/red-orange (more distinct from amber)
        'marker': '^',           # Triangle up
        'size': 38,              # Substantially smaller - secondary
        'edgewidth': 0.4,
        'label': 'P/LP Recall',
        'zorder': 3
    },
    'B/LB Recall': {
        'color': '#4682B4',      # Muted medium blue (steel blue)
        'marker': 's',           # Square
        'size': 38,              # Same as other recalls
        'edgewidth': 0.4,
        'label': 'B/LB Recall',
        'zorder': 3
    },
    'VUS Recall': {
        'color': '#DAA520',      # Muted amber/gold (goldenrod)
        'marker': 'D',           # Diamond
        'size': 38,              # Same as other recalls
        'edgewidth': 0.4,
        'label': 'VUS Recall',
        'zorder': 3
    }
}

metrics = ['Weighted F1 (3-class)', 'P/LP Recall', 'B/LB Recall', 'VUS Recall']

# Create figure
print("\n2. Creating 2x2 dot plot figure...")
fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
axes = axes.flatten()

# Create custom legend handles with clearly visible sizes (independent of plot sizes)
legend_handles = []
legend_marker_sizes = {
    'Weighted F1 (3-class)': 10,  # Larger for primary metric
    'P/LP Recall': 8,             # Smaller for secondary metrics
    'B/LB Recall': 8,
    'VUS Recall': 8
}

for metric in metrics:
    style = metric_styles[metric]
    handle = plt.Line2D([0], [0], marker=style['marker'], color='w',
                       markerfacecolor=style['color'], 
                       markeredgecolor='white',
                       markeredgewidth=0.5,
                       markersize=legend_marker_sizes[metric],
                       label=style['label'], 
                       linestyle='None')
    legend_handles.append(handle)

# Process each cohort
for idx, (cohort_name, cohort_info) in enumerate(cohorts.items()):
    ax = axes[idx]
    
    if cohort_name not in cohort_data:
        ax.text(0.5, 0.5, f'Data not available\nfor {cohort_name}', 
               ha='center', va='center', fontsize=11, color='#999999')
        ax.axis('off')
        continue
    
    df = cohort_data[cohort_name].copy()
    
    # Prepare data
    tools = df['Tool'].tolist()
    n_tools = len(tools)
    
    # Y positions (reversed so highest F1 is at top)
    y_pos = np.arange(n_tools)[::-1]
    
    # Plot each metric with larger markers and white edges
    # Skip NaN values (undefined metrics for classes with zero support)
    for metric in metrics:
        values = df[metric].values
        style = metric_styles[metric]
        
        # Filter out NaN values and corresponding y positions
        valid_mask = ~np.isnan(values)
        valid_values = values[valid_mask]
        valid_y_pos = y_pos[valid_mask]
        
        if len(valid_values) > 0:
            ax.scatter(valid_values, valid_y_pos, 
                      s=style['size'],
                      c=style['color'],
                      marker=style['marker'],
                      alpha=0.95,  # Higher alpha for better visibility
                      edgecolors='white',
                      linewidths=style['edgewidth'],
                      zorder=style['zorder'],
                      label=style['label'])
    
    # Styling
    ax.set_yticks(y_pos)
    ax.set_yticklabels(tools, fontsize=9)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])
    
    # Subtle horizontal grid lines
    ax.grid(axis='x', alpha=0.25, linestyle='--', linewidth=0.5, color='#CCCCCC')
    ax.set_axisbelow(True)
    
    # Subtle horizontal guides at each tool
    for y in y_pos:
        ax.axhline(y=y, color='#F0F0F0', linewidth=0.5, zorder=1)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    
    # Panel label and title
    panel_label = cohort_info['panel_label']
    title = cohort_info['title']
    n = cohort_info['n']
    
    ax.text(-0.12, 1.05, panel_label, transform=ax.transAxes,
           fontsize=16, fontweight='bold', va='top', ha='right')
    ax.text(0.5, 1.03, f'{title} (n={n:,})', transform=ax.transAxes,
           fontsize=10, fontweight='bold', ha='center', va='bottom')
    
    # X-axis label (only bottom row)
    if idx >= 2:
        ax.set_xlabel('Performance Score', fontsize=10, fontweight='bold', color='#333333')
    else:
        ax.set_xlabel('')

print("\n3. Adding shared legend...")

# Adjust layout - optimal spacing with legend closer to panels
plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.subplots_adjust(left=0.10, right=0.98, top=0.89, bottom=0.06, 
                   wspace=0.30, hspace=0.28)

# Add legend at top center with clear separation from panel titles
fig.legend(handles=legend_handles, 
          loc='upper center',
          bbox_to_anchor=(0.5, 0.96),
          ncol=4,
          frameon=False,
          fontsize=11,           # Clear readable font
          columnspacing=3.0,     # Generous space between columns
          handletextpad=0.9,     # Good space between marker and text
          borderpad=1.0,
          markerscale=1.0)       # Use exact custom marker sizes

# Save
print("\n4. Saving figure...")

output_png = output_dir / 'Figure4_Performance_DotPlot.png'
output_pdf = output_dir / 'Figure4_Performance_DotPlot.pdf'

plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, dpi=300, bbox_inches='tight', facecolor='white')

print("\n" + "=" * 80)
print("✓ FIGURE 4 CREATED (DOT PLOT DESIGN)")
print("=" * 80)
print(f"\nOutput files:")
print(f"  - {output_png}")
print(f"  - {output_pdf}")
print(f"\nDesign features:")
print(f"  ✓ Horizontal multi-marker dot plots")
print(f"  ✓ Weighted F1 (3-class) as large primary marker")
print(f"  ✓ P/LP, B/LB, VUS Recall as smaller markers")
print(f"  ✓ Call Rate removed (covered in Figure 2)")
print(f"  ✓ Colorblind-friendly palette (Wong's colors)")
print(f"  ✓ Distinct marker shapes for each metric")
print(f"  ✓ Tools sorted by Weighted F1 (highest at top)")
print(f"  ✓ Subtle horizontal guides")
print(f"  ✓ Clean, minimal design")
print(f"  ✓ Artifact rows removed (analysis, Genebe_Gene)")
print(f"  ✓ HGMD–ClinVar panel titles")
print("=" * 80)

# Print summary statistics
print("\n" + "=" * 80)
print("SUMMARY STATISTICS")
print("=" * 80)

for cohort_name, df in cohort_data.items():
    print(f"\n{cohort_name}:")
    print(f"  Tools: {len(df)}")
    print(f"  Mean Weighted F1: {df['Weighted F1 (3-class)'].mean():.3f}")
    print(f"  Top performer: {df.iloc[0]['Tool']} (F1={df.iloc[0]['Weighted F1 (3-class)']:.3f})")
    print(f"  Mean P/LP Recall: {df['P/LP Recall'].mean():.3f}")
    print(f"  Mean B/LB Recall: {df['B/LB Recall'].mean():.3f}")
    print(f"  Mean VUS Recall: {df['VUS Recall'].mean():.3f}")

print("\n" + "=" * 80)
