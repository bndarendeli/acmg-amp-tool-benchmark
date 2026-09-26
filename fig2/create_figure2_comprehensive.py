#!/usr/bin/env python3
"""
Figure 2: Comprehensive Variant Classification Tool Performance Analysis
Nature Genetics Standard (Arial font, 300 DPI, bold A/B/C panel labels)
All data dynamically calculated from merged_results.tsv - NO HARDCODED VALUES

Optimized layout: 2x2 grid (A left, B top-right, C bottom-right)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Nature Genetics Typography
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 8

print("=" * 100)
print("FIGURE 2: COMPREHENSIVE VARIANT CLASSIFICATION ANALYSIS - OPTIMIZED")
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
# PANEL A DATA: CALL RATE SPECTRUM (FROM CSV)
# ==============================================================================
print("\n2. Loading Panel A: Call Rate Spectrum from CSV...")

# Load corrected call rates from CSV
call_rates_csv = '/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last/results/corrected_metrics/corrected_call_rates.csv'
call_rates_df = pd.read_csv(call_rates_csv)

# Pivot to get tools x datasets
panel_a_df = call_rates_df.pivot(index='Tool', columns='Dataset', values='Call_Rate').reset_index()

# Rename columns to match expected names
panel_a_df.columns.name = None
if 'ClinGen' not in panel_a_df.columns:
    panel_a_df = panel_a_df.rename(columns={'clingen_28012026': 'ClinGen'})
if 'HGMD+ClinVar HL' not in panel_a_df.columns:
    panel_a_df = panel_a_df.rename(columns={'HGMD_Clinvar_HL': 'HGMD+ClinVar HL'})
if 'HGMD+ClinVar Cancer' not in panel_a_df.columns:
    panel_a_df = panel_a_df.rename(columns={'HGMD_Clinvar_Cancer': 'HGMD+ClinVar Cancer'})

# Ensure all expected columns exist
datasets_cols = ['ClinGen', 'FOXL2', 'HGMD+ClinVar HL', 'HGMD+ClinVar Cancer']
for col in datasets_cols:
    if col not in panel_a_df.columns:
        panel_a_df[col] = 0.0

# Calculate mean and sort by Mean (average across all datasets)
panel_a_df['Mean'] = panel_a_df[datasets_cols].mean(axis=1)
panel_a_df = panel_a_df.sort_values(by='Mean', ascending=True).reset_index(drop=True)

print(f"   Tools analyzed: {len(panel_a_df)}")
print(f"   Sample FOXL2 values: {panel_a_df[panel_a_df['Tool'].str.contains('InterVar')][['Tool', 'FOXL2']].to_dict('records')}")

# ==============================================================================
# PANEL B DATA: UNKNOWN VARIANTS BY GROUND TRUTH
# ==============================================================================
print("\n3. Calculating Panel B: Unknown Variants by Ground Truth...")

# Select specific tools for Panel B: Exomiser, InterVar_2018, InterVar_2025, Franklin, DiabloACMG, AutoGVP
# Show both Unknown and NULL/Empty distributions
tools_for_panel_b = ['Exomiser', 'InterVar_2018', 'InterVar_2025', 'Franklin', 'DiabloACMG', 'AutoGVP']

# For each tool in ClinGen, get ground truth distribution of Unknown AND NULL variants
panel_b_data = []

for tool in tools_for_panel_b:
    classification_col = f'{tool}_Classification'
    if classification_col not in clingen_df.columns:
        continue
    
    # Unknown mask
    unknown_mask = (clingen_df[classification_col] == 'Unknown')
    unknown_count = unknown_mask.sum()
    
    # NULL/Empty mask
    null_mask = (clingen_df[classification_col].isna()) | (clingen_df[classification_col] == '')
    null_count = null_mask.sum()
    
    # Total unresolved (Unknown + NULL)
    total_unresolved = unknown_count + null_count
    
    # Ground truth distribution of unknowns
    if unknown_count > 0:
        gt_dist_unknown = clingen_df[unknown_mask]['Ground_Truth_Classification'].value_counts()
        unk_p = gt_dist_unknown.get('Pathogenic', 0)
        unk_lp = gt_dist_unknown.get('Likely_Pathogenic', 0)
        unk_vus = gt_dist_unknown.get('Uncertain_Significance', 0)
        unk_lb = gt_dist_unknown.get('Likely_Benign', 0)
        unk_b = gt_dist_unknown.get('Benign', 0)
    else:
        unk_p = unk_lp = unk_vus = unk_lb = unk_b = 0
    
    # Ground truth distribution of NULLs
    if null_count > 0:
        gt_dist_null = clingen_df[null_mask]['Ground_Truth_Classification'].value_counts()
        null_p = gt_dist_null.get('Pathogenic', 0)
        null_lp = gt_dist_null.get('Likely_Pathogenic', 0)
        null_vus = gt_dist_null.get('Uncertain_Significance', 0)
        null_lb = gt_dist_null.get('Likely_Benign', 0)
        null_b = gt_dist_null.get('Benign', 0)
    else:
        null_p = null_lp = null_vus = null_lb = null_b = 0
    
    # Only include tools with unresolved variants
    if total_unresolved > 0:
        panel_b_data.append({
            'Tool': tool,
            'Total_Unresolved': total_unresolved,
            'Unknown_P': unk_p,
            'Unknown_LP': unk_lp,
            'Unknown_VUS': unk_vus,
            'Unknown_LB': unk_lb,
            'Unknown_B': unk_b,
            'NULL_P': null_p,
            'NULL_LP': null_lp,
            'NULL_VUS': null_vus,
            'NULL_LB': null_lb,
            'NULL_B': null_b
        })

panel_b_df = pd.DataFrame(panel_b_data)
# Sort by total unresolved count descending
panel_b_df = panel_b_df.sort_values('Total_Unresolved', ascending=False).reset_index(drop=True)

print(f"   Tools with unresolved variants: {len(panel_b_df)}")
print(f"   Tools included: {panel_b_df['Tool'].tolist()}")

# Panel C will show Unknown and NULL counts (no separate data calculation needed)

# ==============================================================================
# CREATE FIGURE - OPTIMIZED 2x2 LAYOUT
# ==============================================================================
print("\n5. Creating figure...")

# Figure setup: 2x2 grid, optimized for space
fig = plt.figure(figsize=(16, 10))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.2, 1], height_ratios=[1, 1], 
              hspace=0.25, wspace=0.3)

# ==============================================================================
# PANEL A: CALL RATE SPECTRUM (LEFT, SPANS 2 ROWS)
# ==============================================================================
print("   Creating Panel A: Call Rate Spectrum...")

# Panel A spans both rows on the left
gs_a = gs[:, 0].subgridspec(1, 2, width_ratios=[0.3, 14], wspace=0.01)
ax_a_label = fig.add_subplot(gs_a[0])
ax_a = fig.add_subplot(gs_a[1])

# Panel A label
ax_a_label.text(-5, 1, 'a', fontsize=20, fontweight='bold', 
                ha='center', va='center', fontname='Arial')
ax_a_label.axis('off')

# Plot grouped bars
n_tools = len(panel_a_df)
bar_height = 0.18
spacing = 0.05
y_positions = np.arange(n_tools)

# Dataset colors (matching original)
colors = {
    'ClinGen': '#2E5090',
    'FOXL2': '#C4515C',
    'HGMD+ClinVar HL': '#4DBEAA',
    'HGMD+ClinVar Cancer': '#6B6B6B'
}

# Plot bars for each dataset
for i, dataset in enumerate(datasets_cols):
    offset = (i - 1.5) * (bar_height + spacing)
    values = panel_a_df[dataset].values
    ax_a.barh(y_positions + offset, values, height=bar_height, 
              color=colors[dataset], label=dataset, zorder=3)

# Styling
ax_a.set_yticks(y_positions)
ax_a.set_yticklabels(panel_a_df['Tool'], fontsize=9, fontname='Arial')
ax_a.set_xlabel('Call Rate (%)', fontsize=11, fontweight='bold', fontname='Arial')
ax_a.set_xlim(0, 105)
ax_a.grid(axis='x', alpha=0.2, zorder=1)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# Legend
legend = ax_a.legend(title='Evaluation Cohorts', loc='lower right', 
                     frameon=False, ncol=2, fontsize=9,
                     title_fontproperties={'weight': 'bold', 'size': 10})

# ==============================================================================
# PANEL B: UNKNOWN VARIANTS BY GROUND TRUTH (TOP RIGHT)
# ==============================================================================
print("   Creating Panel B: Unknown Variants by Ground Truth...")

# Panel B in top right
gs_b = gs[0, 1].subgridspec(1, 2, width_ratios=[0.3, 14], wspace=0.01)
ax_b_label = fig.add_subplot(gs_b[0])
ax_b = fig.add_subplot(gs_b[1])

# Panel B label
ax_b_label.text(-5, 1, 'b', fontsize=20, fontweight='bold', 
                ha='center', va='center', fontname='Arial')
ax_b_label.axis('off')

# Vertical stacked bar chart: Ground truth distribution of ALL unresolved variants (Unknown + NULL combined)
n_tools_b = len(panel_b_df)
x_pos_b = np.arange(n_tools_b)

# Nature Communications color-blind friendly red-to-green spectrum
gt_colors = {
    'Pathogenic': '#CA3542',           # Strong red - high pathogenicity
    'Likely_Pathogenic': '#F4A582',    # Salmon/coral - moderate pathogenicity  
    'VUS': '#FFFFBF',                  # Pale yellow - uncertain
    'Likely_Benign': '#92C5DE',        # Light blue-cyan - moderate benign
    'Benign': '#2C7BB6'                # Strong blue - benign
}

# Plot stacked bars - Combined Unknown + NULL (no hatching)
bottom = np.zeros(n_tools_b)

# Plot combined Unknown + NULL for each category
for category in ['Pathogenic', 'Likely_Pathogenic', 'VUS', 'Likely_Benign', 'Benign']:
    # Sum Unknown and NULL for this category
    if category == 'Pathogenic':
        values = panel_b_df['Unknown_P'].values + panel_b_df['NULL_P'].values
    elif category == 'Likely_Pathogenic':
        values = panel_b_df['Unknown_LP'].values + panel_b_df['NULL_LP'].values
    elif category == 'VUS':
        values = panel_b_df['Unknown_VUS'].values + panel_b_df['NULL_VUS'].values
    elif category == 'Likely_Benign':
        values = panel_b_df['Unknown_LB'].values + panel_b_df['NULL_LB'].values
    else:  # Benign
        values = panel_b_df['Unknown_B'].values + panel_b_df['NULL_B'].values
    
    ax_b.bar(x_pos_b, values, bottom=bottom, 
             color=gt_colors[category], width=0.7,
             edgecolor='black', linewidth=0.5)
    bottom += values

# Styling
ax_b.set_xticks(x_pos_b)
ax_b.set_xticklabels(panel_b_df['Tool'], rotation=25, ha='right', fontsize=8, fontname='Arial')
ax_b.set_ylabel('Unresolved Variants (count)', fontsize=10, 
                fontweight='bold', fontname='Arial')
#ax_b.set_title('Ground Truth of Unresolved Variants (Unknown + NULL)', fontsize=10, 
#               fontweight='bold', fontname='Arial', pad=10)
ax_b.grid(axis='y', alpha=0.2, zorder=1)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# Legend: Only classification (P, LP, VUS, LB, B) - removed Unknown/NULL
class_handles = [
    mpatches.Patch(facecolor=gt_colors['Pathogenic'], edgecolor='black', label='P'),
    mpatches.Patch(facecolor=gt_colors['Likely_Pathogenic'], edgecolor='black', label='LP'),
    mpatches.Patch(facecolor=gt_colors['VUS'], edgecolor='black', label='VUS'),
    mpatches.Patch(facecolor=gt_colors['Likely_Benign'], edgecolor='black', label='LB'),
    mpatches.Patch(facecolor=gt_colors['Benign'], edgecolor='black', label='B')
]

legend_b = ax_b.legend(handles=class_handles, loc='upper right', frameon=False, fontsize=8,
                       title_fontsize=8)

# ==============================================================================
# PANEL C: EXPLICIT UNRESOLVED AND NULL/ NO RETURN CLASSIFICATION COUNTS (BOTTOM RIGHT)
# ==============================================================================
print("   Creating Panel C: Explicit Unresolved and NULL/No Return Classification Counts...")

# Panel C in bottom right
gs_c = gs[1, 1].subgridspec(1, 2, width_ratios=[0.3, 14], wspace=0.01)
ax_c_label = fig.add_subplot(gs_c[0])
ax_c = fig.add_subplot(gs_c[1])

# Panel C label
ax_c_label.text(-5, 1, 'c', fontsize=20, fontweight='bold', 
                ha='center', va='center', fontname='Arial')
ax_c_label.axis('off')

# Grouped bar chart: Explicit Unresolved vs NULL/No Return Classification counts
x_pos_c = np.arange(len(panel_b_df))
bar_width = 0.35

# Colors for Explicit Unresolved vs NULL/No Return Classification (matching color scheme)
unknown_color = '#808284'  
null_color = '#D2D7D6'     # Light blue-cyan (matching Likely_Benign)

# Explicit Unresolved bars
ax_c.bar(x_pos_c - bar_width/2, panel_b_df['Total_Unresolved'] - panel_b_df[['NULL_P', 'NULL_LP', 'NULL_VUS', 'NULL_LB', 'NULL_B']].sum(axis=1), 
         bar_width, label='Explicit Unresolved', color=unknown_color, edgecolor='black', linewidth=0.5)

# NULL/No Return Classification bars
ax_c.bar(x_pos_c + bar_width/2, panel_b_df[['NULL_P', 'NULL_LP', 'NULL_VUS', 'NULL_LB', 'NULL_B']].sum(axis=1), 
         bar_width, label='NULL / No returned classification', color=null_color, edgecolor='black', linewidth=0.5)

# Add count labels on bars
for i, row in panel_b_df.iterrows():
    unknown_count = row['Total_Unresolved'] - (row['NULL_P'] + row['NULL_LP'] + row['NULL_VUS'] + row['NULL_LB'] + row['NULL_B'])
    null_count = row['NULL_P'] + row['NULL_LP'] + row['NULL_VUS'] + row['NULL_LB'] + row['NULL_B']
    
    if unknown_count > 0:
        ax_c.text(i - bar_width/2, unknown_count + 50, 
                  str(int(unknown_count)), ha='center', va='bottom', 
                  fontsize=7, fontname='Arial')
    if null_count > 0:
        ax_c.text(i + bar_width/2, null_count + 50, 
                  str(int(null_count)), ha='center', va='bottom', 
                  fontsize=7, fontname='Arial')

# Styling
ax_c.set_xticks(x_pos_c)
ax_c.set_xticklabels(panel_b_df['Tool'], rotation=25, ha='right', fontsize=8, fontname='Arial')
ax_c.set_ylabel('Variant Count', fontsize=10, fontweight='bold', fontname='Arial')
#ax_c.set_title('Unresolved Variant Categories', fontsize=10, fontweight='bold', 
#               fontname='Arial', pad=10)
ax_c.legend(loc='upper right', frameon=False, fontsize=8)
ax_c.grid(axis='y', alpha=0.2, zorder=1)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# ==============================================================================
# SAVE FIGURE
# ==============================================================================
print("\n6. Saving figure...")

# Output to analysis_scripts/figures directory
script_dir = Path(__file__).parent
output_dir = script_dir.parent / 'figures'
output_dir.mkdir(parents=True, exist_ok=True)

output_png = output_dir / 'Figure2_Comprehensive.png'
output_pdf = output_dir / 'Figure2_Comprehensive.pdf'

plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, dpi=300, bbox_inches='tight', facecolor='white')

print("\n" + "=" * 100)
print("✓ FIGURE 2 CREATED - COMPREHENSIVE ANALYSIS")
print("=" * 100)
print(f"\nOutput files:")
print(f"  - {output_png}")
print(f"  - {output_pdf}")
print(f"\nLayout: 2x2 Grid")
print(f"  A (Left, spans 2 rows): Call Rate Spectrum ({len(panel_a_df)} tools, {len(datasets_cols)} cohorts)")
print(f"  B (Top Right): Ground Truth of Unresolved Variants ({len(panel_b_df)} tools)")
print(f"  C (Bottom Right): Explicit Unresolved and NULL/No Return Classification Counts ({len(panel_b_df)} tools)")
print(f"\nAll data dynamically calculated from: {data_file}")
print("=" * 100)
