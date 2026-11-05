#!/usr/bin/env python3
"""
Plot Clustered Pairwise Pearson Correlation Heatmap (Figure 5)

Displays the original 28 x 28 scope: the 25 criteria with nonzero variance
(PP5, BP3 and BP6 omitted) plus MJ, TJ and 3CA, using the original blue palette
and layout.

Statistics are NOT recalculated here. Values are read from the corrected primary
outputs of `calculate_pairwise_pearson_correlations.py`:
  - eligible criterion-criterion pairs: Pearson r, '*' if family A BH-FDR q < 0.05
  - eligible criterion-aggregate and aggregate-aggregate pairs: Pearson r,
    '*' if family B BH-FDR q < 0.05
  - any pair involving a sparse (ineligible) criterion: 'NE' on a neutral gray
    background (not evaluated in the primary correlation analysis)

Primary eligibility (exploratory): nonzero variance and positive Jaccard in >= 5/14
instances; 3/14 and 7/14 are retained as sensitivity analyses in `sensitivity/`.

Hierarchical clustering (average linkage, distance = 1 - Pearson r) uses only the
valid matrix (eligible criteria + aggregates). Sparse criteria are shown
in a separate block in conventional ACMG/AMP order and are not part of the dendrogram.
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
SENSITIVITY_TEXT = '3/14, 4/14 and 7/14'
ALPHA = 0.05
EXPECTED_ZERO_VAR = ['PP5', 'BP3', 'BP6']
EXPECTED_ELIGIBLE = ['PVS1', 'PS1', 'PS3', 'PM1', 'PM2', 'PM4', 'PM5', 'PP2', 'PP3',
                     'BA1', 'BS1', 'BS2', 'BP1', 'BP4', 'BP7']
EXPECTED_SPARSE = ['PS2', 'PS4', 'PM3', 'PM6', 'PP1', 'PP4',
                   'BS3', 'BS4', 'BP2', 'BP5']

# Original Figure 5 palette and mapping (unchanged)
PALETTE = ["#E3F2FD", "#90CAF9", "#0D47A1"]
VMIN, VMAX = 0, 1
NE_COLOR = '#D9D9D9'      # neutral gray, not part of the Pearson colour scale
NE_TEXT_COLOR = '#616161'
DIVIDER_COLOR = '#757575'

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
zero_var = [c for c in ALL_CRITERIA_ORDER if sd[c] == 0]
eligible_rule = [c for c in ALL_CRITERIA_ORDER if n_pos[c] >= MIN_POSITIVE and sd[c] > 0]
displayed_criteria = [c for c in ALL_CRITERIA_ORDER if c not in zero_var]
sparse = [c for c in displayed_criteria if c not in eligible_rule]

assert zero_var == EXPECTED_ZERO_VAR, f"Zero-variance criteria: {zero_var}"
assert eligible_rule == EXPECTED_ELIGIBLE, f"Eligible criteria: {eligible_rule}"
assert sparse == EXPECTED_SPARSE, f"Sparse criteria: {sparse}"
assert list(qc.index[qc['Eligible_primary']]) == eligible_rule, \
    f"QC file disagrees with the {MIN_POSITIVE}/{n_obs} rule"
assert set(long['Analysis']) == {f'primary_min{MIN_POSITIVE}of{n_obs}'}, \
    f"Correlation outputs were not generated with the {MIN_POSITIVE}/{n_obs} rule: {set(long['Analysis'])}"
assert (qc.loc[ALL_CRITERIA_ORDER, 'N_positive_jaccard'] == n_pos[ALL_CRITERIA_ORDER]).all()
assert len(displayed_criteria) == 25
print(f"   Zero variance (omitted): {', '.join(zero_var)}")
print(f"   Eligible ({len(eligible_rule)}): {', '.join(eligible_rule)}")
print(f"   Sparse, displayed as NE ({len(sparse)}): {', '.join(sparse)}")

# ---------------------------------------------------------------------------
# 2. Validate corrected statistical outputs
# ---------------------------------------------------------------------------
print("\n2. Validating corrected statistical outputs...")
valid_vars = eligible_rule + AGGREGATES
assert list(df_r.index) == list(df_r.columns) == valid_vars
assert list(df_fdr.index) == list(df_fdr.columns) == valid_vars
n_elig = len(eligible_rule)
EXPECTED_N = {FAMILY_A: n_elig * (n_elig - 1) // 2, FAMILY_B: n_elig * len(AGGREGATES) + 3}
assert EXPECTED_N == {FAMILY_A: 105, FAMILY_B: 48}, EXPECTED_N
fam_counts = long['FDR_family'].value_counts()
assert all(fam_counts.get(f, 0) == n for f, n in EXPECTED_N.items()), fam_counts
assert set(long['FDR_family']) == {FAMILY_A, FAMILY_B}
tested_vars = set(long['Variable_1']) | set(long['Variable_2'])
assert tested_vars == set(valid_vars), "Ineligible variables present in FDR families"
assert not tested_vars & set(sparse) and not tested_vars & set(zero_var)
for fam, n in EXPECTED_N.items():
    sub = long[long['FDR_family'] == fam]
    assert (sub['FDR_family_size'] == n).all()
    # Independent BH recomputation within the family must reproduce the saved q-values
    q_check = multipletests(sub['P_value'], method='fdr_bh')[1]
    assert np.allclose(q_check, sub['FDR_q_value'], rtol=0, atol=1e-15), f"BH mismatch in {fam}"
    assert ((q_check < ALPHA) == sub['FDR_significant'].values).all()
fam_lookup = {}
for row in long.itertuples():
    assert df_r.loc[row.Variable_1, row.Variable_2] == row.Pearson_r == df_r.loc[row.Variable_2, row.Variable_1]
    assert df_fdr.loc[row.Variable_1, row.Variable_2] == row.FDR_q_value == df_fdr.loc[row.Variable_2, row.Variable_1]
    assert bool(row.FDR_significant) == (row.FDR_q_value < ALPHA)
    fam_lookup[frozenset((row.Variable_1, row.Variable_2))] = row.FDR_family
assert np.allclose(np.diag(df_r.values), 1.0)

# r values must be the full-sample (n = 14) Pearson r of the unchanged master table
for row in long.itertuples():
    r_check = np.corrcoef(master[row.Variable_1], master[row.Variable_2])[0, 1]
    assert np.isclose(r_check, row.Pearson_r, rtol=0, atol=1e-12), (row.Variable_1, row.Variable_2)
print(f"   ✓ Family A: {fam_counts[FAMILY_A]} tests; family B: {fam_counts[FAMILY_B]} tests")
print("   ✓ No sparse or zero-variance criterion in either FDR family")
print("   ✓ r matches full-sample Pearson r from the master table (n = 14)")

# The original palette maps only [0, 1]; stop rather than misrepresent negatives
off = ~np.eye(len(valid_vars), dtype=bool)
min_r = float(df_r.values[off].min())
if min_r < 0:
    neg = long.loc[long['Pearson_r'] < 0, ['Variable_1', 'Variable_2', 'Pearson_r']]
    raise ValueError("Negative eligible correlations found; the original 0-1 blue palette cannot "
                     f"represent them without clipping. Review the design:\n{neg.to_string(index=False)}")
print(f"   ✓ No negative eligible correlations (minimum off-diagonal r = {min_r:.3f})")

# ---------------------------------------------------------------------------
# 3. Clustering on the valid 16-variable matrix only
# ---------------------------------------------------------------------------
print("\n3. Hierarchical clustering (eligible variables only)...")
print("   Distance: 1 - Pearson r; linkage: average")
distance_matrix = 1 - df_r.values.astype(float)
np.fill_diagonal(distance_matrix, 0)
linkage_matrix = linkage(squareform(distance_matrix, checks=False), method='average')
clustered_valid = [valid_vars[i] for i in dendrogram(linkage_matrix, no_plot=True)['leaves']]
VARIABLE_ORDER = clustered_valid + sparse
n_valid, n_vars = len(clustered_valid), len(VARIABLE_ORDER)
assert n_vars == 28 and len(set(VARIABLE_ORDER)) == 28
assert not set(VARIABLE_ORDER) & set(zero_var)
display_labels = [AGGREGATE_LABELS.get(v, v) for v in VARIABLE_ORDER]

# ---------------------------------------------------------------------------
# 4. Display-only matrix (NE for any pair involving a sparse criterion)
# ---------------------------------------------------------------------------
display = pd.DataFrame('NE', index=VARIABLE_ORDER, columns=VARIABLE_ORDER, dtype=object)
r_plot = np.full((n_vars, n_vars), np.nan)
for i, v1 in enumerate(VARIABLE_ORDER):
    for j, v2 in enumerate(VARIABLE_ORDER):
        if v1 in sparse or v2 in sparse:
            continue
        r_plot[i, j] = df_r.loc[v1, v2]
        if i == j:
            display.iloc[i, j] = ''
        else:
            q = df_fdr.loc[v1, v2]
            display.iloc[i, j] = f"{r_plot[i, j]:.2f}" + ('*' if q < ALPHA else '')
for s in sparse:
    display.loc[s, s] = ''
assert (display.values == display.values.T).all()
assert np.array_equal(np.isnan(r_plot), np.isnan(r_plot.T))
assert np.allclose(r_plot[~np.isnan(r_plot)], r_plot.T[~np.isnan(r_plot.T)])
display.index = display.columns = display_labels
display.to_csv(output_dir / 'figure5_display_matrix.csv')

# ---------------------------------------------------------------------------
# 5. Figure (original layout)
# ---------------------------------------------------------------------------
print("\n4. Creating clustered heatmap with dendrogram...")
fig = plt.figure(figsize=(14, 13))
gs = fig.add_gridspec(1, 2, width_ratios=[0.6, 5], wspace=0.08)

ax_dend = fig.add_subplot(gs[0, 0])
linkage_plot = linkage_matrix.copy()
linkage_plot[:, 2] = np.power(linkage_plot[:, 2], 0.35)  # display-only branch stretch (as original)
dendrogram(linkage_plot, orientation='left', ax=ax_dend, color_threshold=0,
           above_threshold_color='black', no_labels=True)
for line in ax_dend.get_lines():
    line.set_linewidth(1.3)
for collection in ax_dend.collections:
    collection.set_linewidth(1.3)
# Leaf k is drawn at y = 10k + 5; invert and extend so the 16 leaves align with
# heatmap rows 0..15 and the sparse block (rows 16..27) has no branches.
ax_dend.set_ylim(10 * n_vars, 0)
ax_dend.set_xticks([])
ax_dend.set_yticks([])
for spine in ax_dend.spines.values():
    spine.set_visible(False)
ax_dend.text(0.5, 10 * (n_valid + len(sparse) / 2), 'Sparse criteria\n(not clustered)',
             transform=ax_dend.get_yaxis_transform(), ha='center', va='center',
             fontsize=8, style='italic', color=NE_TEXT_COLOR)

ax_heat = fig.add_subplot(gs[0, 1])
cmap = LinearSegmentedColormap.from_list("pearson_custom", PALETTE, N=256)
norm = Normalize(vmin=VMIN, vmax=VMAX)
# Background layer: neutral gray for NE cells (outside the Pearson colour scale)
ax_heat.imshow(np.ones((n_vars, n_vars, 3)) * np.array(matplotlib.colors.to_rgb(NE_COLOR)),
               aspect='auto', interpolation='nearest')
im = ax_heat.imshow(np.ma.masked_invalid(r_plot), cmap=cmap, norm=norm,
                    aspect='auto', interpolation='nearest')

cbar = plt.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04, aspect=20)
cbar.set_label('Pearson r', rotation=270, labelpad=20, fontsize=11)
cbar.ax.tick_params(labelsize=9)

ax_heat.set_xticks(np.arange(n_vars))
ax_heat.set_yticks(np.arange(n_vars))
ax_heat.set_xticklabels(display_labels, rotation=90, ha='center', fontsize=8)
ax_heat.set_yticklabels(display_labels, fontsize=8)
ax_heat.tick_params(axis='y', which='major', pad=2)
for tick in ax_heat.get_xticklabels() + ax_heat.get_yticklabels():
    if tick.get_text() in sparse:
        tick.set_color(NE_TEXT_COLOR)
        tick.set_style('italic')

n_markers = {FAMILY_A: 0, FAMILY_B: 0}
n_ne = 0
for i, v1 in enumerate(VARIABLE_ORDER):
    for j, v2 in enumerate(VARIABLE_ORDER):
        if i == j:
            continue
        text = display.iloc[i, j]
        if text == 'NE':
            n_ne += 1
            ax_heat.text(j, i, 'NE', ha='center', va='center', fontsize=6, color=NE_TEXT_COLOR)
            continue
        r = r_plot[i, j]
        if text.endswith('*') and i < j:
            n_markers[fam_lookup[frozenset((v1, v2))]] += 1
        ax_heat.text(j, i, text, ha='center', va='center', fontsize=6,
                     color='white' if abs(r) > 0.5 else 'black', fontweight='normal')

for pos in [n_valid - 0.5]:
    ax_heat.axhline(pos, color=DIVIDER_COLOR, linewidth=1.0)
    ax_heat.axvline(pos, color=DIVIDER_COLOR, linewidth=1.0)

plt.tight_layout(rect=[0, 0.085, 1, 1])

pdf_path = output_dir / 'Figure5_clustered_pairwise_pearson_heatmap.pdf'
png_path = output_dir / 'Figure5_clustered_pairwise_pearson_heatmap.png'
plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
plt.savefig(png_path, dpi=300, bbox_inches='tight')
plt.close()

# ---------------------------------------------------------------------------
# 6. Post-plot validation
# ---------------------------------------------------------------------------
n_sparse_pairs = n_vars * (n_vars - 1) - n_valid * (n_valid - 1)
assert n_ne == n_sparse_pairs, (n_ne, n_sparse_pairs)
for fam in [FAMILY_A, FAMILY_B]:
    expected = int(long.loc[long['FDR_family'] == fam, 'FDR_significant'].sum())
    assert n_markers[fam] == expected, (fam, n_markers[fam], expected)
pos_of = {v: k for k, v in enumerate(VARIABLE_ORDER)}
for fam, sub in long.groupby('FDR_family'):
    q_check = multipletests(sub['P_value'], method='fdr_bh')[1]
    for row, q in zip(sub.itertuples(), q_check):
        i, j = pos_of[row.Variable_1], pos_of[row.Variable_2]
        for cell in (display.iloc[i, j], display.iloc[j, i]):
            assert cell == f"{row.Pearson_r:.2f}" + ('*' if q < ALPHA else ''), (row.Variable_1, row.Variable_2)
for s in sparse:
    k = pos_of[s]
    assert set(display.iloc[k, :]) <= {'NE', ''} and set(display.iloc[:, k]) <= {'NE', ''}

print(f"   Dimensions: {n_vars} x {n_vars} ({n_valid} clustered + {len(sparse)} sparse)")
print(f"   Significance markers (unique pairs): family A {n_markers[FAMILY_A]}/{EXPECTED_N[FAMILY_A]}, "
      f"family B {n_markers[FAMILY_B]}/{EXPECTED_N[FAMILY_B]}")
print("   ✓ Every plotted marker agrees with BH-FDR recomputed within its family")
print(f"   NE cells (off-diagonal): {n_ne}")
print(f"   ✓ PDF: {pdf_path}")
print(f"   ✓ PNG: {png_path}")
print(f"   ✓ Display-only matrix: {output_dir / 'figure5_display_matrix.csv'}")
print("\nDisplay order:")
for i, lab in enumerate(display_labels, 1):
    print(f"   {i:2d}. {lab}{'  (sparse, NE)' if VARIABLE_ORDER[i - 1] in sparse else ''}")
print("=" * 80)
