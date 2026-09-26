# ACMG/AMP variant interpretation benchmark: data and analysis archive

Supporting analysis files for a benchmark of 12 ACMG/AMP interpretation tools (14 benchmark instances) across ClinGen, FOXL2, HGMD–ClinVar hearing-loss and cancer-predisposition cohorts. The supplied manuscript draft describes 15,433 variants.

This repository contains merged predictions and reference annotations, 52 matching raw tool-output candidates, derived tables, analysis scripts and figure exports. Raw outputs are losslessly gzip-compressed; original and stored SHA-256 hashes are recorded. The current analysis layout covers Figures 2–6; historical manifests and reproduction helpers still refer to the earlier layout.

**Before publication:** the archived ClinGen Figure 3 Jaccard table includes 257 variants with missing reference criteria, whereas the current script and method notes exclude them. This affects 109 of 392 tool–criterion count pairs. Values have not been silently corrected; see the [input verification report](docs/INPUT_VERIFICATION.md).

## Start here

- [Data catalog](docs/DATA_CATALOG.md): every CSV, its row count and column names.
- [Input catalog](docs/INPUT_CATALOG.md): processed TSV schemas and links to the selected raw outputs.
- [Figure and supplementary-table guide](docs/FIGURE_GUIDE.md): navigate the manuscript-related outputs.
- [Reproduction guide](docs/REPRODUCIBILITY.md): available inputs, executable commands and missing dependencies.
- [Tool-result release plan](docs/RESULTS_RELEASE.md): which prediction files and metadata are needed for reproducibility.
- [Verified input selection](docs/INPUT_VERIFICATION.md): matches, remaining interpretation issues and execution evidence.
- [Archive notes](docs/ARCHIVE_NOTES.md): draft metadata, historical inconsistencies and release scope.
- [Validation results](docs/VALIDATION.md): completed checks and runtime limitations.
- [Push instructions](docs/GITHUB_PUSH.md).

## Repository layout

```text
data/
  processed/   Frozen merged predictions and reference classes/criteria
  tool_results/<tool>/<cohort>/   Matching raw-output candidates (.gz)
analysis/
  fig2/       Classification availability, call rates and tool parsers
  fig3/       Criterion counts and criterion-level Jaccard analysis
  fig4/       Classification performance and corrected metric tables
  fig5/       Pairwise Pearson correlations: data/, scripts/, figure/
  fig6/       Evidence concordance, capability adjustment: outputs/, scripts/, figures/
docs/         Catalog, provenance, figure guide and reproduction notes
metadata/     Input checksums, source selection and audit results
tools/        Archive validation and isolated figure reproduction
requirements.txt
```

Each analysis directory contains its own figure exports. The former shared figures directory has been removed. Figure 5 now contains the clustered Pearson correlation analysis. Evidence-concordance outputs are under Figure 6; some retain historical `fig4` or `figure5` names. Supporting documentation may still reference the earlier layout.

## Figures and outputs

| Figure | Analysis | Main PDF |
| --- | --- | --- |
| 2 | Classification availability and call rates | [Figure 2](analysis/fig2/figures/Figure2_Comprehensive.pdf) |
| 3 | Criterion counts and criterion-level Jaccard | [ClinGen Figure 3](analysis/fig3/figures/Figure3_ClinGen_Counts_Jaccard.pdf) |
| 4 | Reconstructed three-class accuracy and class-specific F1 | [Figure 4](analysis/fig4/figures/Figure4_Performance_DotPlot.pdf) |
| 5 | Clustered pairwise Pearson correlations | [Figure 5](analysis/fig5/figure/Figure5_clustered_pairwise_pearson_heatmap.pdf) |
| 6 | Evidence concordance and capability adjustment | [Figure 6 preview](analysis/fig6/figures/Figure6_five_panel_preview_v3.pdf) |

Figure 5 uses 14 tool instances as observations and 31 numeric variables: 28 criterion-level Jaccard measures, mean criterion Jaccard, total criterion Jaccard and three-class accuracy. The plot omits three zero-variance criteria (PP5, BP3 and BP6). The master table and Pearson correlation, p-value and FDR tables are in [Figure 5 data](analysis/fig5/data/).

[Figure 6 outputs](analysis/fig6/outputs/) include the variant–tool master table, ClinGen and FOXL2 capability-adjusted summaries, total criterion Jaccard values and tool capability tables. Supplementary figures are stored alongside the relevant main figures; their filenames retain historical numbering.

## Input scope

The main processed TSV contains 15,433 rows: 15,145 labeled main-cohort rows and 288 unlabeled rows. The separate FOXL2 TSV contains 288 rows. Use explicit cohort filters rather than blindly concatenating the two files. See the [data overview](data/README.md).

## Setup and reproduction status

```bash
python -m pip install -r requirements.txt
python -m pip install scipy statsmodels
```

Figure 5 imports SciPy and statsmodels, which are not explicitly listed in the current requirements file. Dependencies are unpinned. Historical runs used Python 3.12.11 under Ubuntu WSL2; see the [runtime versions](metadata/runtime_versions.json) and [WSL notes](docs/LOCAL_WSL.md).

**Reproduction helpers still target the earlier Figures 2–5 layout.** The current Figures 2–6 package has the following path limitations:

- `tools/reproduce_figures.py --figure 4` refers to the absent `create_figure4_dotplot.py`; the current script is `analysis/fig4/create_figure4.py`.
- Figure 3 plotting expects the master table under `analysis/fig5/outputs/`; it is now under `analysis/fig6/outputs/`.
- `tools/reproduce_from_inputs.py` references earlier Figure 5 scripts that are no longer present in that directory.
- Figure 5 scripts retain input paths that do not match the current `data/` and `scripts/` split. Figure 6 scripts also retain historical or working-directory-dependent paths.

These paths and required inputs must be reconciled before rerunning the full analysis. The historical [reproduction guide](docs/REPRODUCIBILITY.md) does not establish successful execution of the reorganized package.

The existing raw-output reparsing helper is:

```bash
python tools/rebuild_tool_predictions.py --output reproduced/raw-rerun
```

This reconstructs predictions from saved outputs, not the original tool/database runs. See the input verification report for historical comparisons and limitations. Reproduction outputs belong under ignored `reproduced/`.

## Validate the archive

Requires Python 3.9+; no third-party packages are needed for validation.

```bash
python tools/validate_repository.py
```

The checker examines recorded hashes, Python syntax, table structure and file sizes. Its manifests describe the earlier archive, and it explicitly rejects `analysis/fig6/`. It therefore needs updating before a clean validation result can be expected for the current layout. The linked validation report records checks from 2026-09-12, not verification of this reorganization.

The historical ClinGen Figure 3 missing-reference discrepancy remains documented above; reorganizing files does not resolve it. Validation does not establish scientific conclusions or exact agreement with the manuscript.

## Citation and reuse

The supplied PDF still has placeholder title and author metadata. A publication DOI, final title, author list and repository URL were not supplied, so no bibliographic metadata has been invented. Cite the final associated article when available.

No code or data license was selected on the authors' behalf. Third-party source-data terms are not replaced by this archive. The files include outputs and derived records for HGMD–ClinVar benchmark cohorts; packaging does not establish their redistribution terms. Original tool databases and the manuscript PDF are not included.

