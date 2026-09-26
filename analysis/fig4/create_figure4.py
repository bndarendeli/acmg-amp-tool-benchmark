#!/usr/bin/env python3
"""
Figure 4: Classification Performance - DOT PLOT DESIGN
2x2 Multi-Panel Layout with Horizontal Multi-Marker Dot Plots

UPDATED: 2026-09-21
Uses reconstructed metrics with 3-class Accuracy (NOT Weighted F1)

Design:
- Accuracy (3-class): Large filled circle (primary metric, purple)
- P/LP F1: Triangle (orange/red)
- VUS F1: Diamond (gold)
- B/LB F1: Square (blue)
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
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
print("FIGURE 4: CLASSIFICATION PERFORMANCE (RECONSTRUCTED METRICS)")
print("=" * 80)

# Paths
base_dir = Path(__file__).parent
data_dir = base_dir / 'data' / 'reconstructed'
output_dir = base_dir / 'figures'
output_dir.mkdir(parents=True, exist_ok=True)

# Load reconstructed metrics
print("\n1. Loading reconstructed metrics...")

cohorts = {
    'ClinGen': {
        'n': 11409,
        'file': 'ClinGen_reconstructed_metrics.csv',
        'panel_label': 'a',
        'title': 'ClinGen Master Cohort',
        'expected_tools': 14
    },
    'FOXL2': {
        'n': 288,
        'file': 'FOXL2_reconstructed_metrics.csv',
        'panel_label': 'b',
        'title': 'FOXL2 Disease-Specific Cohort',
        'expected_tools': 12
    },
    'HGMD–ClinVar HL': {
        'n': 1948,
        'file': 'HGMD_ClinVar_HL_reconstructed_metrics.csv',
        'panel_label': 'c',
        'title': 'HGMD–ClinVar Hearing Loss',
        'expected_tools': 14
    },
    'HGMD–ClinVar Cancer': {
        'n': 1788,
        'file': 'HGMD_ClinVar_Cancer_reconstructed_metrics.csv',
        'panel_label': 'd',
        'title': 'HGMD–ClinVar Cancer Predisposition',
        'expected_tools': 14
    }
}

# Metric definitions
metrics = {
    'accuracy_3class': {
        'label': 'Accuracy (3-class)',
        'color': '#4B0082',  # Deep indigo/purple
        'marker': 'o',
        'size': 85
    },
    'P_LP_f1': {
        'label': 'P/LP F1',
        'color': '#D94F30',  # Vermillion/red-orange
        'marker': '^',
        'size': 45
    },
    'VUS_f1': {
        'label': 'VUS F1',
        'color': '#DAA520',  # Amber/gold
        'marker': 'D',
        'size': 45
    },
    'B_LB_f1': {
        'label': 'B/LB F1',
        'color': '#4682B4',  # Steel blue
        'marker': 's',
        'size': 45
    }
}

# Load and validate data
cohort_data = {}
all_data = []

for cohort_name, cohort_info in cohorts.items():
    file_path = data_dir / cohort_info['file']
    
    print(f"\n  Loading {cohort_name}...")
    df = pd.read_csv(file_path)
    
    # Validation 1: Tool count
    n_tools = len(df)
    expected = cohort_info['expected_tools']
    print(f"    Tools: {n_tools} (expected {expected})")
    assert n_tools == expected, f"Tool count mismatch for {cohort_name}"
    
    # Validation 2: Cohort size
    assert (df['cohort_n'] == cohort_info['n']).all(), f"Cohort size mismatch for {cohort_name}"
    
    # Validation 3: Accuracy calculation
    accuracy_check = np.isclose(df['accuracy_3class'], df['n_correct_3class'] / df['cohort_n'])
    assert accuracy_check.all(), f"Accuracy calculation error in {cohort_name}"
    
    # Validation 4: Value ranges
    for metric_col in ['accuracy_3class', 'P_LP_f1', 'VUS_f1', 'B_LB_f1']:
        if metric_col in df.columns:
            values = df[metric_col].dropna()
            assert ((values >= 0) & (values <= 1)).all(), f"Invalid {metric_col} values in {cohort_name}"
    
    # Validation 5: FOXL2 B/LB check
    if cohort_name == 'FOXL2':
        assert df['B_LB_f1'].isna().all(), "FOXL2 should have NA for B_LB_f1"
        print(f"    ✓ FOXL2 B_LB_f1 all NA (no B/LB reference variants)")
    
    # Validation 6: No Weighted F1 column
    assert 'Weighted F1 (3-class)' not in df.columns, "Weighted F1 column found - should not exist"
    
    # Sort by accuracy (descending)
    df = df.sort_values('accuracy_3class', ascending=False).reset_index(drop=True)
    
    # Assert highest accuracy is in first row
    assert df.iloc[0]['accuracy_3class'] == df['accuracy_3class'].max(), \
        f"Sorting error in {cohort_name}: first row should have max accuracy"
    
    cohort_data[cohort_name] = df
    
    # Add to combined data
    df_copy = df.copy()
    df_copy['cohort'] = cohort_name
    all_data.append(df_copy)
    
    print(f"    ✓ Validation passed")

print("\n✓ All validations passed")

# Create combined CSV
print("\n2. Creating combined classification-performance CSV...")
df_combined = pd.concat(all_data, ignore_index=True)

# Reorder columns
column_order = [
    'cohort', 'tool', 'cohort_n', 'n_valid_predictions', 'n_correct_3class',
    'accuracy_3class', 'P_LP_f1', 'VUS_f1', 'B_LB_f1',
    'P_LP_precision', 'P_LP_recall', 'P_LP_support',
    'VUS_precision', 'VUS_recall', 'VUS_support',
    'B_LB_precision', 'B_LB_recall', 'B_LB_support'
]
df_combined = df_combined[column_order]

combined_csv_path = output_dir / 'Figure4_combined_classification_performance.csv'
df_combined.to_csv(combined_csv_path, index=False)
print(f"  ✓ Saved: {combined_csv_path}")

# Create figure
print("\n3. Creating Figure 4...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()


# Collect legend handles from first panel
legend_handles = []
legend_labels = []

for idx, (cohort_name, cohort_info) in enumerate(cohorts.items()):
    ax = axes[idx]
    df = cohort_data[cohort_name]
    
    n_tools = len(df)
    # CRITICAL FIX: Invert y_positions so highest accuracy (row 0) appears at top
    y_positions = np.arange(n_tools)[::-1]
    
    # Verify first row has max accuracy and will be plotted at top
    assert df.iloc[0]['accuracy_3class'] == df['accuracy_3class'].max(), \
        f"First row should have max accuracy in {cohort_name}"
    
    # Plot each metric
    for metric_name, metric_info in metrics.items():
        # Skip B/LB for FOXL2
        if cohort_name == 'FOXL2' and metric_name == 'B_LB_f1':
            continue
        
        values = df[metric_name].values
        
        scatter = ax.scatter(values, y_positions,
                            color=metric_info['color'],
                            marker=metric_info['marker'],
                            s=metric_info['size'],
                            alpha=0.9,
                            edgecolor='white',
                            linewidth=0.5,
                            zorder=5 if metric_name == 'accuracy_3class' else 3,
                            label=metric_info['label'])
        
        # Collect legend handles from first panel
        if idx == 0 and cohort_name == 'ClinGen':
            legend_handles.append(scatter)
            legend_labels.append(metric_info['label'])
    
    # Y-axis: tool names
    ax.set_yticks(y_positions)
    ax.set_yticklabels(df['tool'], fontsize=8)
    ax.tick_params(axis='y', length=3, pad=4)
    
    # X-axis
    if idx >= 2:
        ax.set_xlabel('Performance Score', fontsize=10, fontweight='bold')
    else:
        ax.set_xlabel('')
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(axis='x', labelsize=8)
    
    # Grid
    ax.grid(True, axis='x', alpha=0.3, linewidth=0.5, linestyle='-', color='gray')
    ax.set_axisbelow(True)
    
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_linewidth(0.8)
    
    # Panel label (bold, outside)
    ax.text( -0.15, 1.05, cohort_info['panel_label'],
            transform=ax.transAxes,
            fontsize=14, fontweight='bold',
            va='top', ha='left')
    
    # Title
    ax.set_title(
        f"{cohort_info['title']} (n={cohort_info['n']:,})",
        fontsize=11,
        pad=10,
        fontweight='bold',
        loc='center'
    )
plt.tight_layout(rect=[0, 0, 1, 0.95])

# Add shared legend at top center
fig.legend(legend_handles, legend_labels,
          loc='upper center',
          ncol=4,
          fontsize=10,
          frameon=False,
          bbox_to_anchor=(0.5, 1.00),
          columnspacing=2.0,
          handletextpad=0.5)

# Save figure
print("\n4. Saving figure...")

png_path = output_dir / 'Figure4_Performance_DotPlot.png'
pdf_path = output_dir / 'Figure4_Performance_DotPlot.pdf'

plt.savefig(png_path, dpi=300, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')

print(f"  ✓ PNG: {png_path}")
print(f"  ✓ PDF: {pdf_path}")

plt.close()

# Summary
print("\n" + "=" * 80)
print("FIGURE 4 GENERATION COMPLETE")
print("=" * 80)

print("\nValidation Results:")
print(f"  ✓ ClinGen: {len(cohort_data['ClinGen'])} tools")
print(f"  ✓ FOXL2: {len(cohort_data['FOXL2'])} tools")
print(f"  ✓ HGMD–ClinVar HL: {len(cohort_data['HGMD–ClinVar HL'])} tools")
print(f"  ✓ HGMD–ClinVar Cancer: {len(cohort_data['HGMD–ClinVar Cancer'])} tools")
print(f"  ✓ All accuracy values = n_correct / cohort_n")
print(f"  ✓ All values in range [0, 1]")
print(f"  ✓ FOXL2 B_LB_f1 all NA")
print(f"  ✓ No Weighted F1 column used")

print("\nOutput Files:")
print(f"  ✓ {png_path}")
print(f"  ✓ {pdf_path}")
print(f"  ✓ {combined_csv_path}")

print("\n" + "=" * 80)
