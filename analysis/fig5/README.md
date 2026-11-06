# Figure 5

Clustered Pearson correlation heatmap of criterion-level Jaccard concordance profiles
and aggregate performance metrics across the 14 ClinGen tool instances.

Output: `Figure5_clustered_pairwise_pearson_heatmap.pdf` (and `.png`, 300 dpi).

## Reproduce

The four Figure 5 steps were verified in Ubuntu WSL using the existing
`acmg-benchmarking` Conda environment with Python 3.9.23. Check the interpreter before
running: `conda run -n acmg-benchmarking python -c "import sys; print(sys.executable); print(sys.version)"`.
Use that environment for all four steps; avoid mixing Windows Python with WSL Python.

`environment.yml` records the Python and direct package versions used in the successful
Figure 5 check. To create a separate Figure 5 environment without changing an existing
one, run `conda env create -f analysis/fig5/environment.yml`, then
`conda activate acmg-benchmarking-fig5`. This file covers Figure 5 only; the complete
repository pipeline also needs the packages in the root `requirements.txt`.

From the repository root, activate the chosen environment and run in order:

```bash
python analysis/fig5/build_correlation_master_table.py
python analysis/fig5/calculate_pairwise_pearson_correlations.py
python analysis/fig5/compare_criterion_assignment_counts_with_jaccard.py   # descriptive QC only
python analysis/fig5/plot_clustered_pairwise_pearson_heatmap.py
```

Each script writes into this directory and can be called from any working directory.
`python tools/run_all.py` runs the same four steps, after the upstream Figure 3, 4 and 6
steps, in a separate copy under `reproduced/`; add `--from-raw` to start from the bundled
tool outputs.

| Step | Script | Reads | Writes |
|---|---|---|---|
| 1 | `build_correlation_master_table.py` | `fig3/jaccard_per_criterion_clingen_28012026.csv`, `fig6/outputs/fig6_capability_adjusted_summary.csv` (MJ), `fig6/outputs/clingen_total_criterion_jaccard.csv` (TJ), `fig4/data/reconstructed/ClinGen_reconstructed_metrics.csv` (3CA) | `correlation_master_table.csv` (14 tool instances x 32 columns) |
| 2 | `calculate_pairwise_pearson_correlations.py` | master table (read only; checked against the step-1 sources) | QC, primary matrices, long table, summary, `sensitivity/` |
| 3 | `compare_criterion_assignment_counts_with_jaccard.py` | `data/processed/merged_results_noHGMD.tsv`, `fig6/outputs/fig6_master_variant_tool.csv`, Figure 3 Jaccard table | `criterion_assignment_vs_jaccard_*.csv` (not used by the figure) |
| 4 | `plot_clustered_pairwise_pearson_heatmap.py` | master table, QC file, primary matrices, long table | Figure 5 PDF and PNG |

## Analysis

- Observations: 14 tool instances; all 14 values, including zeros, are used for every
  correlation. Source Jaccard values and aggregate metrics are not modified or recomputed.
- Eligibility: positive criterion-level Jaccard in >= 5/14 instances and nonzero
  variance (`criterion_correlation_qc.csv`). This is an exploratory threshold, not an
  established or prospectively prespecified cutoff. Results at 3/14, 4/14 and 7/14 are
  evaluated in sensitivity analyses.
- Eligible criteria (15): PVS1, PS1, PS3, PM1, PM2, PM4, PM5, PP2, PP3, BA1, BS1, BS2,
  BP1, BP4, BP7. Not eligible: PS2, PS4, PM3, PM6, PP1, PP4, BS3, BS4, BP2, BP5
  (positive Jaccard in fewer than 5 instances) and PP5, BP3, BP6 (zero variance).
- Pearson r with two-sided p; Benjamini–Hochberg FDR in two separate families:
  A = 105 criterion–criterion pairs; B = 45 criterion–aggregate + 3 aggregate–aggregate
  pairs (48 tests). MJ = mean criterion Jaccard, TJ = total criterion Jaccard,
  3CA = three-class accuracy.

## Figure

18 x 18 matrix: the 15 eligible criteria plus MJ, TJ and 3CA; ineligible criteria are
not shown. Cells give Pearson r (two decimals); `*` marks BH-FDR q < 0.05 within the
pair's testing family. Rows and columns are ordered by average-linkage hierarchical
clustering on 1 − r (display order only; dendrogram branch heights are stretched for
readability). Colours use the blue scale from 0 to 1; the script stops if any
correlation is negative, because this scale cannot show it. Before plotting, the script
recomputes BH within each family and checks every annotation against the saved tables.

## Outputs

Primary (used by the figure): `criterion_correlation_qc.csv`, `pearson_r_matrix.csv`,
`pearson_pvalue_matrix.csv`, `pearson_fdr_matrix.csv` (family-specific q),
`pearson_fdr_family_matrix.csv`, `pearson_pairwise_long.csv` (with leave-one-out and
tool-dependence robustness columns), `correlation_analysis_summary.md`.

Descriptive only: `criterion_assignment_vs_jaccard_sparsity.csv` and
`criterion_assignment_vs_jaccard_summary.csv` (zero Jaccard vs. never assigned).

Sensitivity (`sensitivity/`): `sensitivity_threshold_*.csv` (>= 3, 4, 5, 7 of 14),
`sensitivity_leave_one_out_*_min{5,7}of14.csv` and `sensitivity_tool_family_*_min{5,7}of14.csv`
(leave-one-out and one-representative-per-tool-family analyses for the 5/14 primary and
the 7/14 eligibility sets), `tool_instance_families.csv`,
`tool_instance_profile_similarity.csv`, and `sensitivity_nonzero_variance_pairwise_long.csv`
(378-test, single-family sensitivity analysis, for comparison only).
