# Figure 5

Selected output: `Figure5_clustered_pairwise_pearson_heatmap.pdf` (and `.png`).

Run, in order (from any working directory):

```bash
python analysis/fig5/build_correlation_master_table.py
python analysis/fig5/calculate_pairwise_pearson_correlations.py
python analysis/fig5/compare_criterion_assignment_counts_with_jaccard.py
python analysis/fig5/plot_clustered_pairwise_pearson_heatmap.py
```

`correlation_master_table.csv` (14 tool instances x 32 columns) is read only by the
correlation script; values are checked against the Figure 3/4/6 source tables.

## Analysis

- Criterion sparsity QC: `criterion_correlation_qc.csv`.
- Primary eligibility: positive criterion-level Jaccard in >= 5/14 instances and nonzero
  variance; all 14 observations, including zeros, are kept. This is an exploratory
  threshold, not an established or prospectively prespecified cutoff. Results at 3/14,
  4/14 and 7/14 are retained as sensitivity analyses.
- 15 criteria are eligible. BH-FDR families: A = 105 criterion–criterion pairs among
  eligible criteria; B = 45 criterion–aggregate pairs plus the three aggregate–aggregate
  pairs (MJ, TJ, 3CA), 48 tests.

## Figure display

The heatmap keeps the original 28 x 28 scope (25 criteria with nonzero variance plus
MJ, TJ, 3CA; PP5, BP3 and BP6 omitted) and the original blue palette (0-1). Pairs
involving one of the 10 sparse criteria are shown as `NE` (not evaluated in the primary
correlation analysis) on a gray background and are in neither FDR family. Clustering
(average linkage, 1 - r) uses only the 18 eligible variables; sparse criteria form a
separate block in ACMG/AMP order. The plot script stops if any eligible r is negative,
because the 0-1 palette cannot show it. `figure5_display_matrix.csv` is display-only.

## Outputs

Primary (used by the figure): `pearson_r_matrix.csv`, `pearson_pvalue_matrix.csv`,
`pearson_fdr_matrix.csv` (family-specific q), `pearson_fdr_family_matrix.csv`,
`pearson_pairwise_long.csv` (with leave-one-out and dependence robustness columns),
`correlation_analysis_summary.md`.

Descriptive only: `criterion_assignment_vs_jaccard_sparsity.csv` and
`criterion_assignment_vs_jaccard_summary.csv` (zero Jaccard vs. never assigned).

Sensitivity (`sensitivity/`): `sensitivity_threshold_*.csv` (>= 3, 4, 5, 7 of 14),
`sensitivity_leave_one_out_*_min{5,7}of14.csv` and `sensitivity_tool_family_*_min{5,7}of14.csv`
(leave-one-out and one-representative-per-tool-family analyses for the 5/14 primary and
the retained 7/14 eligibility sets),
`tool_instance_families.csv`, `tool_instance_profile_similarity.csv`, and
`legacy_zero_variance_only_pairwise_long.csv` (reproduction of the original 378-test,
single-family analysis, for comparison only).
