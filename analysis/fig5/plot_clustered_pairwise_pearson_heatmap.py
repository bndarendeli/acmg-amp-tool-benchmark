#!/usr/bin/env python3
"""
Plot Clustered Pairwise Pearson Correlation Heatmap (Figure 5)

Displays the primary correlation matrix produced by
`calculate_pairwise_pearson_correlations.py`: the eligible ACMG/AMP criteria
(positive criterion-level Jaccard in >= 5/14 tool instances and nonzero variance)
plus MJ, TJ and 3CA. Criteria that do not meet the eligibility rule are not shown.

Statistics are NOT recalculated here; values are read from the saved outputs and
validated against them:
  - cell values: Pearson r across all 14 tool instances (two decimals)
  - '*': BH-FDR q < 0.05 within the corresponding testing family
      family A = criterion-criterion pairs
      family B = criterion-aggregate and aggregate-aggregate pairs
Hierarchical clustering (average linkage, distance = 1 - Pearson r) affects the
display order only. Colours use the Figure 5 blue palette (r from 0 to 1).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests

print("=" * 80)
print("PLOTTING CLUSTERED PAIRWISE PEARSON CORRELATION HEATMAP (FIGURE 5)")
print("=" * 80)

input_dir = Path(__file__).resolve().parent
output_dir = input_dir
paths = {
    'master': input_dir / 'correlation_master_table.csv',
    'qc': input_dir / 'criterion_correlation_qc.csv',
    'r': input_dir / 'pearson_r_matrix.csv',
    'fdr': input_dir / 'pearson_fdr_matrix.csv',
    'long': input_dir / 'pearson_pairwise_long.csv',
}
for path in paths.values():
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path} "
                                "(run calculate_pairwise_pearson_correlations.py first)")

ALL_CRITERIA_ORDER = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]
AGGREGATES = ['Mean_Criterion_Jaccard', 'Total_Criterion_Jaccard', 'Accuracy_3class']
AGGREGATE_LABELS = {'Mean_Criterion_Jaccard': 'MJ', 'Total_Criterion_Jaccard': 'TJ',
                    'Accuracy_3class': '3CA'}
FAMILY_A = 'A_criterion_criterion'
FAMILY_B = 'B_criterion_aggregate_and_aggregate_aggregate'
MIN_POSITIVE = 5   # exploratory primary threshold; 3/14, 4/14 and 7/14 are sensitivity analyses
ALPHA = 0.05
EXPECTED_ELIGIBLE = ['PVS1', 'PS1', 'PS3', 'PM1', 'PM2', 'PM4', 'PM5', 'PP2', 'PP3',
                     'BA1', 'BS1', 'BS2', 'BP1', 'BP4', 'BP7']

# Figure 5 palette and mapping
PALETTE = ["#E3F2FD", "#90CAF9", "#0D47A1"]
VMIN, VMAX = 0, 1

# ---------------------------------------------------------------------------
# 1. Load inputs and verify eligibility programmatically
# ---------------------------------------------------------------------------
print("\n1. Loading inputs and verifying eligibility...")
master = pd.read_csv(paths['master'])
qc = pd.read_csv(paths['qc']).set_index('Criterion')
df_r = pd.read_csv(paths['r'], index_col=0)
df_fdr = pd.read_csv(paths['fdr'], index_col=0)
long = pd.read_csv(paths['long'])
assert master.shape == (14, 32), f"Master table shape {master.shape}"
n_obs = len(master)

n_pos = (master[ALL_CRITERIA_ORDER] > 0).sum()
sd = master[ALL_CRITERIA_ORDER].std()
eligible = [c for c in ALL_CRITERIA_ORDER if n_pos[c] >= MIN_POSITIVE and sd[c] > 0]
assert eligible == EXPECTED_ELIGIBLE, f"Eligible criteria: {eligible}"
assert list(qc.index[qc['Eligible_primary']]) == eligible, \
    f"QC file disagrees with the {MIN_POSITIVE}/{n_obs} rule"
assert set(long['Analysis']) == {f'primary_min{MIN_POSITIVE}of{n_obs}'}, \
    f"Correlation outputs were not generated with the {MIN_POSITIVE}/{n_obs} rule: {set(long['Analysis'])}"
print(f"   Eligible criteria ({len(eligible)}): {', '.join(eligible)}")

# ---------------------------------------------------------------------------
# 2. Validate the saved statistical outputs
# ---------------------------------------------------------------------------
print("\n2. Validating saved statistical outputs...")
variables = eligible + AGGREGATES
n_vars = len(variables)
assert n_vars == 18
assert list(df_r.index) == list(df_r.columns) == variables
assert list(df_fdr.index) == list(df_fdr.columns) == variables
assert np.array_equal(df_r.values, df_r.values.T), "r matrix not symmetric"
off = ~np.eye(n_vars, dtype=bool)
assert np.array_equal(df_fdr.values[off], df_fdr.values.T[off]), "FDR matrix not symmetric"
assert np.allclose(np.diag(df_r.values), 1.0)

n_elig = len(eligible)
EXPECTED_N = {FAMILY_A: n_elig * (n_elig - 1) // 2, FAMILY_B: n_elig * len(AGGREGATES) + 3}
assert EXPECTED_N == {FAMILY_A: 105, FAMILY_B: 48}, EXPECTED_N
assert set(long['FDR_family']) == set(EXPECTED_N)
assert set(long['Variable_1']) | set(long['Variable_2']) == set(variables)
fam_lookup = {}
for fam, n in EXPECTED_N.items():
    sub = long[long['FDR_family'] == fam]
    assert len(sub) == n and (sub['FDR_family_size'] == n).all()
    # Independent BH recomputation within the family must reproduce the saved q-values
    q_check = multipletests(sub['P_value'], method='fdr_bh')[1]
    assert np.allclose(q_check, sub['FDR_q_value'], rtol=0, atol=1e-15), f"BH mismatch in {fam}"
    assert ((q_check < ALPHA) == sub['FDR_significant'].values).all()
for row in long.itertuples():
    assert df_r.loc[row.Variable_1, row.Variable_2] == row.Pearson_r == df_r.loc[row.Variable_2, row.Variable_1]
    assert df_fdr.loc[row.Variable_1, row.Variable_2] == row.FDR_q_value
    r_check = np.corrcoef(master[row.Variable_1], master[row.Variable_2])[0, 1]
    assert np.isclose(r_check, row.Pearson_r, rtol=0, atol=1e-12), (row.Variable_1, row.Variable_2)
    fam_lookup[frozenset((row.Variable_1, row.Variable_2))] = row.FDR_family
print(f"   ✓ Family A: {EXPECTED_N[FAMILY_A]} tests; family B: {EXPECTED_N[FAMILY_B]} tests; "
      "BH q-values reproduced")
print(f"   ✓ r equals full-sample Pearson r from the master table (n = {n_obs})")

# The palette maps only [0, 1]; stop rather than misrepresent negatives
min_r = float(df_r.values[off].min())
if min_r < 0:
    neg = long.loc[long['Pearson_r'] < 0, ['Variable_1', 'Variable_2', 'Pearson_r']]
    raise ValueError("Negative correlations found; the 0-1 blue palette cannot represent them "
                     f"without clipping. Review the design:\n{neg.to_string(index=False)}")
print(f"   ✓ No negative correlations (minimum off-diagonal r = {min_r:.3f})")

# ---------------------------------------------------------------------------
# 3. Hierarchical clustering (display order only)
# ---------------------------------------------------------------------------
print("\n3. Hierarchical clustering (average linkage, distance = 1 - Pearson r)...")
distance_matrix = 1 - df_r.values.astype(float)
np.fill_diagonal(distance_matrix, 0)
linkage_matrix = linkage(squareform(distance_matrix, checks=False), method='average')
order = [variables[i] for i in dendrogram(linkage_matrix, no_plot=True)['leaves']]
df_r_c = df_r.loc[order, order]
df_fdr_c = df_fdr.loc[order, order]
labels = [AGGREGATE_LABELS.get(v, v) for v in order]
assert np.allclose(np.sort(df_r.values.ravel()), np.sort(df_r_c.values.ravel()))

# ---------------------------------------------------------------------------
# 4. Figure
# ---------------------------------------------------------------------------
print("\n4. Creating clustered heatmap with dendrogram...")
fig = plt.figure(figsize=(11, 10))
gs = fig.add_gridspec(1, 2, width_ratios=[0.6, 5], wspace=0.12)

ax_dend = fig.add_subplot(gs[0, 0])
linkage_plot = linkage_matrix.copy()
linkage_plot[:, 2] = np.power(linkage_plot[:, 2], 0.35)  # display-only branch stretch
dendrogram(linkage_plot, orientation='left', ax=ax_dend, color_threshold=0,
           above_threshold_color='black', no_labels=True)
for line in ax_dend.get_lines():
    line.set_linewidth(1.3)
for collection in ax_dend.collections:
    collection.set_linewidth(1.3)
ax_dend.set_ylim(10 * n_vars, 0)  # leaf k at y = 10k + 5, aligned with heatmap row k
ax_dend.set_xticks([])
ax_dend.set_yticks([])
for spine in ax_dend.spines.values():
    spine.set_visible(False)

ax_heat = fig.add_subplot(gs[0, 1])
cmap = LinearSegmentedColormap.from_list("pearson_custom", PALETTE, N=256)
r_values = df_r_c.values.astype(float)
im = ax_heat.imshow(r_values, cmap=cmap, norm=Normalize(vmin=VMIN, vmax=VMAX),
                    aspect='auto', interpolation='nearest')

cbar = plt.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04, aspect=20)
cbar.set_label('Pearson r', rotation=270, labelpad=20, fontsize=12)
cbar.ax.tick_params(labelsize=10)

ax_heat.set_xticks(np.arange(n_vars))
ax_heat.set_yticks(np.arange(n_vars))
ax_heat.set_xticklabels(labels, rotation=90, ha='center', fontsize=10)
ax_heat.set_yticklabels(labels, fontsize=10)
ax_heat.tick_params(axis='y', which='major', pad=2)

plotted = {}
n_markers = {FAMILY_A: 0, FAMILY_B: 0}
for i, v1 in enumerate(order):
    for j, v2 in enumerate(order):
        if i == j:
            continue
        r = r_values[i, j]
        text = f'{0.0 if abs(r) < 0.005 else r:.2f}'
        if df_fdr_c.loc[v1, v2] < ALPHA:
            text += '*'
            if i < j:
                n_markers[fam_lookup[frozenset((v1, v2))]] += 1
        plotted[(v1, v2)] = text
        ax_heat.text(j, i, text, ha='center', va='center', fontsize=7.5,
                     color='white' if abs(r) > 0.5 else 'black')

fig.subplots_adjust(left=0.03, right=0.97, top=0.98, bottom=0.08)

pdf_path = output_dir / 'Figure5_clustered_pairwise_pearson_heatmap.pdf'
png_path = output_dir / 'Figure5_clustered_pairwise_pearson_heatmap.png'
fig.savefig(pdf_path, dpi=300, bbox_inches='tight')
fig.savefig(png_path, dpi=300, bbox_inches='tight')
plt.close(fig)

# ---------------------------------------------------------------------------
# 5. Post-plot validation: every annotation matches the saved statistics
# ---------------------------------------------------------------------------
assert len(plotted) == n_vars * (n_vars - 1)
for fam in EXPECTED_N:
    expected = int(long.loc[long['FDR_family'] == fam, 'FDR_significant'].sum())
    assert n_markers[fam] == expected, (fam, n_markers[fam], expected)
for row in long.itertuples():
    exp = f"{row.Pearson_r:.2f}" + ('*' if row.FDR_q_value < ALPHA else '')
    assert plotted[(row.Variable_1, row.Variable_2)] == plotted[(row.Variable_2, row.Variable_1)] == exp

print(f"   Dimensions: {n_vars} x {n_vars} ({n_elig} criteria + {len(AGGREGATES)} aggregates)")
print(f"   Significance markers (unique pairs): family A {n_markers[FAMILY_A]}/{EXPECTED_N[FAMILY_A]}, "
      f"family B {n_markers[FAMILY_B]}/{EXPECTED_N[FAMILY_B]}")
print("   ✓ Every annotation matches pearson_pairwise_long.csv (r and BH-FDR marker)")
print(f"   ✓ PDF: {pdf_path}")
print(f"   ✓ PNG: {png_path}")
print("\nClustered display order:")
for i, lab in enumerate(labels, 1):
    print(f"   {i:2d}. {lab}")
print("=" * 80)
