#!/usr/bin/env python3
"""
Plot Clustered Pairwise Pearson Correlation Heatmap

Supplementary visualization showing the same 28×28 Pearson correlation matrix
with hierarchical clustering to reveal correlation structure.

Uses EXISTING Pearson correlations and FDR values - does NOT recalculate.
Clustering only affects display order, not statistics.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

print("=" * 80)
print("PLOTTING CLUSTERED PAIRWISE PEARSON CORRELATION HEATMAP")
print("=" * 80)

# Define paths
input_dir = Path(__file__).parent
output_dir = input_dir

r_matrix_path = input_dir / 'pearson_r_matrix.csv'
fdr_matrix_path = input_dir / 'pearson_fdr_matrix.csv'

# Verify files exist
for path in [r_matrix_path, fdr_matrix_path]:
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")

print(f"\n1. Loading EXISTING Pearson correlation matrices...")
print(f"   r: {r_matrix_path}")
print(f"   FDR: {fdr_matrix_path}")

# Load matrices
df_r = pd.read_csv(r_matrix_path, index_col=0)
df_fdr = pd.read_csv(fdr_matrix_path, index_col=0)

print(f"   Loaded: {df_r.shape[0]} × {df_r.shape[1]}")

# Define variable order (ALL 31)
ALL_CRITERIA_ORDER = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]

AGGREGATES = [
    'Mean_Criterion_Jaccard',
    'Total_Criterion_Jaccard',
    'Accuracy_3class'
]

# Zero-variance criteria to omit
ZERO_VAR_CRITERIA = ['PP5', 'BP3', 'BP6']

# Variables to display (25 criteria + 3 aggregates = 28)
CRITERIA_DISPLAY = [c for c in ALL_CRITERIA_ORDER if c not in ZERO_VAR_CRITERIA]
VARIABLE_ORDER = CRITERIA_DISPLAY + AGGREGATES

# Display labels for aggregates
AGGREGATE_LABELS = {
    'Mean_Criterion_Jaccard': 'MJ',
    'Total_Criterion_Jaccard': 'TJ',
    'Accuracy_3class': '3CA'
}

# Create display labels
display_labels_dict = {}
for var in VARIABLE_ORDER:
    if var in AGGREGATE_LABELS:
        display_labels_dict[var] = AGGREGATE_LABELS[var]
    else:
        display_labels_dict[var] = var

print(f"\n2. Data preparation...")
print(f"   Total variables in full matrix: 31")
print(f"   Zero-variance variables (omitted): {len(ZERO_VAR_CRITERIA)}")
print(f"   Variables to display: {len(VARIABLE_ORDER)}")
print(f"     - Criteria: {len(CRITERIA_DISPLAY)}")
print(f"     - Aggregates: {len(AGGREGATES)}")

# Subset matrices to displayed variables
df_r = df_r.loc[VARIABLE_ORDER, VARIABLE_ORDER]
df_fdr = df_fdr.loc[VARIABLE_ORDER, VARIABLE_ORDER]

# Validation
print(f"\n3. Validation...")
assert df_r.shape == (28, 28), f"Expected 28×28, got {df_r.shape}"
assert df_fdr.shape == (28, 28), f"Expected 28×28, got {df_fdr.shape}"
print(f"   ✓ Exactly 28 variables")

# Check PP5, BP3, BP6 are absent
for var in ZERO_VAR_CRITERIA:
    assert var not in df_r.index, f"{var} should be excluded"
print(f"   ✓ PP5, BP3, BP6 are absent")

# Check diagonal is 1.0
diagonal_vals = np.diag(df_r.values)
assert np.allclose(diagonal_vals, 1.0, rtol=1e-5), "Diagonal should be 1.0"
print(f"   ✓ Diagonal is 1.0")

# Check symmetry
assert np.allclose(df_r.values, df_r.values.T, rtol=1e-10), "Matrix should be symmetric"
print(f"   ✓ Matrix is symmetric")

# Perform hierarchical clustering
print(f"\n4. Performing hierarchical clustering...")
print(f"   Distance metric: 1 - Pearson r")
print(f"   Linkage method: average")

# Convert correlation to distance
# For valid correlations, distance = 1 - r
# This ensures perfect positive correlation (r=1) has distance 0
r_array = df_r.values.astype(float)

# Create distance matrix
distance_matrix = 1 - r_array

# Ensure diagonal is exactly 0
np.fill_diagonal(distance_matrix, 0)

# Convert to condensed distance matrix (upper triangle)
condensed_dist = squareform(distance_matrix, checks=False)

# Perform hierarchical clustering
linkage_matrix = linkage(condensed_dist, method='average')

# Get dendrogram to extract order
dend = dendrogram(linkage_matrix, no_plot=True)
clustered_indices = dend['leaves']

# Reorder variables based on clustering
clustered_var_order = [VARIABLE_ORDER[i] for i in clustered_indices]

print(f"   ✓ Clustering complete")

# Reorder matrices
df_r_clustered = df_r.loc[clustered_var_order, clustered_var_order]
df_fdr_clustered = df_fdr.loc[clustered_var_order, clustered_var_order]

# Validation after clustering
print(f"\n5. Post-clustering validation...")
assert df_r_clustered.shape == (28, 28), "Shape changed after clustering"
print(f"   ✓ Shape preserved (28×28)")

# Check row order == column order
assert list(df_r_clustered.index) == list(df_r_clustered.columns), \
    "Row and column order must match"
print(f"   ✓ Clustered row order == clustered column order")

# Check symmetry preserved
assert np.allclose(df_r_clustered.values, df_r_clustered.values.T, rtol=1e-10), \
    "Symmetry lost after clustering"
print(f"   ✓ Matrix remains symmetric after reordering")

# Check diagonal still 1.0
diagonal_clustered = np.diag(df_r_clustered.values)
assert np.allclose(diagonal_clustered, 1.0, rtol=1e-5), "Diagonal changed"
print(f"   ✓ Diagonal remains 1.0")

# Check values unchanged (just reordered)
original_values_sorted = np.sort(df_r.values.flatten())
clustered_values_sorted = np.sort(df_r_clustered.values.flatten())
assert np.allclose(original_values_sorted, clustered_values_sorted, rtol=1e-10), \
    "Values changed during clustering"
print(f"   ✓ All r values identical to original (only reordered)")

# Create display labels in clustered order
clustered_display_labels = [display_labels_dict[var] for var in clustered_var_order]

print(f"\n6. Creating clustered heatmap with dendrogram...")

# Set up figure with dendrogram and heatmap
fig = plt.figure(figsize=(14, 13))

# Create grid: dendrogram on left, heatmap on right
gs = fig.add_gridspec(1, 2, width_ratios=[0.6, 5], wspace=0.08)

# Panel a: Dendrogram
ax_dend = fig.add_subplot(gs[0, 0])

# Keep the original linkage matrix for clustering/order.
# Create a copy only for dendrogram visualization.
linkage_plot = linkage_matrix.copy()

# Square-root transform expands very small linkage distances
# while preserving their ordering.
linkage_plot[:, 2] = np.power(linkage_plot[:, 2], 0.35)

# Plot dendrogram
dend_plot = dendrogram(
    linkage_plot, 
    orientation='left',
    ax=ax_dend,
    color_threshold=0,
    above_threshold_color='black',
    no_labels=True
)

# Make dendrogram branches more visible
for line in ax_dend.get_lines():
    line.set_linewidth(1.3)

for collection in ax_dend.collections:
    collection.set_linewidth(1.3)
    
ax_dend.set_xticks([])
ax_dend.set_yticks([])
for spine in ax_dend.spines.values():
    spine.set_visible(False)
ax_dend.set_ylabel('')
ax_dend.set_xlabel('')

# Panel b: Heatmap
ax_heat = fig.add_subplot(gs[0, 1])

# Create masked array for NAs (shouldn't occur)
r_array_clustered = df_r_clustered.values.astype(float)
r_masked = np.ma.masked_where(np.isnan(r_array_clustered), r_array_clustered)

# Use custom blue-centered diverging colormap
cmap = LinearSegmentedColormap.from_list(
    "pearson_custom",
    [
        "#A82020",  # -1
        "#FFDADA",  # -0.5
        "#FFF2F2",  #  0
        "#A2CB8B",  # +0.5
        "#165823"   # +1
    ],
    N=256
)

# Create heatmap
im = ax_heat.imshow(r_masked, cmap=cmap, aspect='auto', 
                    vmin=-1, vmax=1, interpolation='nearest')

# Add colorbar
cbar = plt.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04)
cbar.set_label('Pearson r', rotation=270, labelpad=20, fontsize=11)
cbar.ax.tick_params(labelsize=9)

# Set ticks and labels
ax_heat.set_xticks(np.arange(len(clustered_display_labels)))
ax_heat.set_yticks(np.arange(len(clustered_display_labels)))
ax_heat.set_xticklabels(clustered_display_labels, rotation=90, ha='center', fontsize=8)
ax_heat.set_yticklabels(clustered_display_labels, fontsize=8)
ax_heat.tick_params(axis='y', which='major', pad=2)

# Add Pearson r values and EXISTING FDR significance markers
print(f"   Adding Pearson r values and EXISTING FDR significance markers...")
n_markers = 0

for i in range(len(clustered_var_order)):
    for j in range(len(clustered_var_order)):
        # Skip diagonal (self-correlations)
        if i == j:
            continue
        
        var1 = clustered_var_order[i]
        var2 = clustered_var_order[j]
        r = df_r_clustered.loc[var1, var2]
        fdr_val = df_fdr_clustered.loc[var1, var2]
        
        # Skip if NA
        if pd.isna(r):
            continue
        
        # Avoid displaying "-0.00"
        r_display = 0.0 if abs(r) < 0.005 else r
        
        # Display Pearson r with two decimals
        text = f'{r_display:.2f}'
        
        # Add marker if EXISTING global FDR significant
        if not pd.isna(fdr_val) and fdr_val < 0.05:
            text += '*'
            n_markers += 1
        
        # Choose text color based on background
        text_color = 'white' if abs(r) > 0.5 else 'black'
        
        ax_heat.text(j, i, text, ha='center', va='center',
                    fontsize=6, color=text_color, fontweight='normal')

print(f"   Significance markers added: {n_markers}")

# Add footnotes
fig.text(0.5, 0.02, 
         '* FDR q < 0.05 (Benjamini–Hochberg correction, global pairwise analysis)\n' +
         'PP5, BP3, and BP6 were omitted because of zero variance.',
         ha='center', fontsize=9, style='italic')

# Adjust layout
plt.tight_layout(rect=[0, 0.04, 1, 1])

# Save outputs
print(f"\n7. Saving outputs...")

pdf_path = output_dir / 'clustered_pairwise_pearson_heatmap.pdf'
png_path = output_dir / 'clustered_pairwise_pearson_heatmap.png'

plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
plt.savefig(png_path, dpi=300, bbox_inches='tight')

print(f"   ✓ PDF: {pdf_path}")
print(f"   ✓ PNG: {png_path}")

plt.close()

# Print clustered variable order
print("\n" + "=" * 80)
print("FINAL CLUSTERED VARIABLE ORDER")
print("=" * 80)
print("\nVariables ordered by hierarchical clustering (average linkage, distance = 1 - r):")
for i, var in enumerate(clustered_var_order, 1):
    display_label = display_labels_dict[var]
    if var != display_label:
        print(f"{i:2d}. {var:<30} ({display_label})")
    else:
        print(f"{i:2d}. {var}")

# Summary
print("\n" + "=" * 80)
print("CLUSTERED PAIRWISE PEARSON HEATMAP COMPLETE")
print("=" * 80)

print(f"\nDimensions: {len(clustered_var_order)} × {len(clustered_var_order)}")
print(f"  - Criteria: {len(CRITERIA_DISPLAY)}")
print(f"  - Aggregates: {len(AGGREGATES)}")
print(f"  - Omitted zero-variance: {', '.join(ZERO_VAR_CRITERIA)}")

print(f"\nClustering:")
print(f"  - Method: Hierarchical (average linkage)")
print(f"  - Distance: 1 - Pearson r")
print(f"  - Effect: Display order only (statistics unchanged)")

print(f"\nSignificance markers: {n_markers} (from EXISTING global FDR analysis)")

print(f"\nOutput files:")
print(f"  ✓ {pdf_path}")
print(f"  ✓ {png_path}")

print(f"\nVisualization notes:")
print(f"  - Symmetric 28×28 matrix with dendrogram")
print(f"  - Same clustered order on both axes")
print(f"  - Correlation scale: -1 to +1, centered at 0")
print(f"  - Custom blue-centered diverging colormap")
print(f"  - Pearson r values in cells (2 decimals)")
print(f"  - Uses EXISTING global FDR q-values (NOT recalculated)")

print("\n" + "=" * 80)
print("No existing files were modified.")
print("Clustering changed only display order, not statistics.")
print("=" * 80)
