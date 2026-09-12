#!/usr/bin/env python3
"""
Supplementary Figure: InterVar 2018 vs 2025 vs Ground Truth Comparison
Grouped bar chart showing classification distributions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Nature Genetics Typography
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10

print("=" * 100)
print("SUPPLEMENTARY FIGURE: INTERVAR 2018 vs 2025 vs GROUND TRUTH")
print("=" * 100)

# ==============================================================================
# LOAD DATA
# ==============================================================================
print("\n1. Loading data...")
data_file = '/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last/data/merged_results.tsv'
df = pd.read_csv(data_file, sep='\t', low_memory=False)
clingen_df = df[df['Dataset'] == 'clingen_28012026'].copy()
total_clingen = len(clingen_df)
print(f"   ClinGen variants: {total_clingen:,}")

# ==============================================================================
# CALCULATE INTERVAR COMPARISON DATA
# ==============================================================================
print("\n2. Calculating InterVar comparison...")

iv2018 = clingen_df['InterVar_2018_Classification'].copy()
iv2025 = clingen_df['InterVar_2025_Classification'].copy()
ground_truth = clingen_df['Ground_Truth_Classification'].copy()

# Classification order for display
class_order = ['Pathogenic', 'Likely_Pathogenic', 'Uncertain_Significance', 
               'Likely_Benign', 'Benign']

# Count distributions
dist_2018 = iv2018.value_counts()
dist_2025 = iv2025.value_counts()
dist_gt = ground_truth.value_counts()

# Create dataframe for plotting
panel_data = []
for cls in class_order:
    count_2018 = dist_2018.get(cls, 0)
    count_2025 = dist_2025.get(cls, 0)
    count_gt = dist_gt.get(cls, 0)
    
    pct_2018 = (count_2018 / total_clingen) * 100
    pct_2025 = (count_2025 / total_clingen) * 100
    pct_gt = (count_gt / total_clingen) * 100
    
    panel_data.append({
        'Classification': cls,
        'InterVar_2018': pct_2018,
        'InterVar_2025': pct_2025,
        'Ground_Truth': pct_gt,
        'Count_2018': count_2018,
        'Count_2025': count_2025,
        'Count_GT': count_gt
    })

# Add Unknown separately (only for InterVar, not ground truth)
unknown_2018 = dist_2018.get('Unknown', 0)
unknown_2025 = dist_2025.get('Unknown', 0)
panel_data.append({
    'Classification': 'Unknown',
    'InterVar_2018': (unknown_2018 / total_clingen) * 100,
    'InterVar_2025': (unknown_2025 / total_clingen) * 100,
    'Ground_Truth': 0,  # Ground truth has no Unknown
    'Count_2018': unknown_2018,
    'Count_2025': unknown_2025,
    'Count_GT': 0
})

panel_df = pd.DataFrame(panel_data)

print(f"   Classification categories: {len(panel_df)}")

# ==============================================================================
# CREATE FIGURE
# ==============================================================================
print("\n3. Creating figure...")

fig, ax = plt.figure(figsize=(10, 6)), plt.gca()

# Bar positions
n_classes = len(panel_df)
x_pos = np.arange(n_classes)
bar_width = 0.25

# Short labels for x-axis
short_labels = {
    'Pathogenic': 'P',
    'Likely_Pathogenic': 'LP',
    'Uncertain_Significance': 'VUS',
    'Likely_Benign': 'LB',
    'Benign': 'B',
    'Unknown': 'Unk'
}

# Different color scheme for supplementary figure (purple/magenta tones)
class_colors = {
    'Pathogenic': '#D53E4F',           # Red-magenta
    'Likely_Pathogenic': '#FC8D59',    # Orange-coral
    'Uncertain_Significance': '#FEE08B',  # Yellow
    'Likely_Benign': '#99D594',        # Light green
    'Benign': '#3288BD',               # Blue
    'Unknown': '#BABABA'               # Light gray
}

# Plot grouped bars (3 bars per classification)
for i, row in panel_df.iterrows():
    cls = row['Classification']
    color = class_colors[cls]
    
    # InterVar 2018 (left) - lighter shade
    ax.bar(i - bar_width, row['InterVar_2018'], bar_width, 
           color=color, alpha=0.4, edgecolor='black', linewidth=0.8)
    
    # InterVar 2025 (middle) - medium shade
    ax.bar(i, row['InterVar_2025'], bar_width, 
           color=color, alpha=0.7, edgecolor='black', linewidth=0.8)
    
    # Ground Truth (right) - full color, only if not Unknown
    if cls != 'Unknown':
        ax.bar(i + bar_width, row['Ground_Truth'], bar_width, 
               color=color, alpha=1.0, edgecolor='black', linewidth=1.2)

# Add count labels on top of bars
for i, row in panel_df.iterrows():
    # 2018 count
    if row['Count_2018'] > 0:
        ax.text(i - bar_width, row['InterVar_2018'] + 0.8, 
                str(int(row['Count_2018'])), ha='center', va='bottom', 
                fontsize=8, fontname='Arial')
    # 2025 count
    if row['Count_2025'] > 0:
        ax.text(i, row['InterVar_2025'] + 0.8, 
                str(int(row['Count_2025'])), ha='center', va='bottom', 
                fontsize=8, fontname='Arial')
    # Ground truth count
    if row['Count_GT'] > 0:
        ax.text(i + bar_width, row['Ground_Truth'] + 0.8, 
                str(int(row['Count_GT'])), ha='center', va='bottom', 
                fontsize=8, fontname='Arial', fontweight='bold')

# Styling
ax.set_xticks(x_pos)
ax.set_xticklabels([short_labels[cls] for cls in panel_df['Classification']], 
                    fontsize=11, fontname='Arial', fontweight='bold')
ax.set_ylabel('Percentage of Variants (%)', fontsize=12, 
              fontweight='bold', fontname='Arial')
#ax.set_title('InterVar 2018 vs 2025 vs Ground Truth Comparison', 
#             fontsize=14, fontweight='bold', fontname='Arial', pad=15)
ax.set_ylim(0, max(panel_df[['InterVar_2018', 'InterVar_2025', 'Ground_Truth']].max()) * 1.15)
ax.grid(axis='y', alpha=0.3, zorder=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend with distinct colors
legend_handles = [
    mpatches.Patch(facecolor='#9467BD', alpha=0.4, edgecolor='black', linewidth=0.8, 
                   label='InterVar 2018'),
    mpatches.Patch(facecolor='#9467BD', alpha=0.7, edgecolor='black', linewidth=0.8, 
                   label='InterVar 2025'),
    mpatches.Patch(facecolor='#9467BD', alpha=1.0, edgecolor='black', linewidth=1.2, 
                   label='Ground Truth')
]
ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1, 1.02), frameon=True, fontsize=10,
          fancybox=True, shadow=True)

plt.tight_layout()

# ==============================================================================
# SAVE FIGURE
# ==============================================================================
print("\n4. Saving figure...")

# Output to analysis_scripts/figures directory
script_dir = Path(__file__).parent
output_dir = script_dir.parent / 'figures'
output_dir.mkdir(parents=True, exist_ok=True)

output_png = output_dir / 'Supp_Fig_InterVar_Comparison.png'
output_pdf = output_dir / 'Supp_Fig_InterVar_Comparison.pdf'

plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, dpi=300, bbox_inches='tight', facecolor='white')

print("\n" + "=" * 100)
print("✓ SUPPLEMENTARY FIGURE CREATED - INTERVAR COMPARISON")
print("=" * 100)
print(f"\nOutput files:")
print(f"  - {output_png}")
print(f"  - {output_pdf}")
print(f"\nData: ClinGen dataset ({total_clingen:,} variants)")
print(f"Classifications: {len(panel_df)} categories")
print("=" * 100)
