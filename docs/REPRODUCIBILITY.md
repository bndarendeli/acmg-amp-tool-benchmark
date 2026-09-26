# Reproducibility

The Figures 2–5 pipeline ran in Ubuntu WSL2, Conda `base`, Python 3.12.11. Installed analysis package versions are recorded in `metadata/runtime_versions.json`; they are not the original ACMG tool environments.

## Validate the package

```bash
python tools/validate_repository.py
```

Checks archived/input hashes, Python syntax and archived CSV structure. The input manifest separately records original hashes of gzip-compressed outputs.

## Reparse raw outputs

```bash
python tools/rebuild_tool_predictions.py --output reproduced/raw-rerun
```

Decompresses selected files temporarily, runs bundled parsers and matches hg19/hg38 keys against reference mappings in the processed TSVs. It retains the first duplicate raw key, matching historical code. Five-tier classifications and canonical criterion sets are compared, treating missing/Unknown as no-calls. The eight FOXL2 VIP-HL encoding differences remain documented. This reconstructs predictions from saved outputs, not the original tool/database runs.

## Recalculate analyses and regenerate figures

```bash
python tools/reproduce_from_inputs.py --output reproduced/analysis-rerun
```

Defaults read `data/processed/merged_results.tsv` and `foxl2_merged_results.tsv`. Optional `--merged` and `--foxl2` flags select other inputs. The runner copies the archive, redirects hardcoded paths/stale parser imports in that copy, recalculates Figure 4 metrics and Figure 5 master/capability/supplementary tables, and runs Figures 2–5 and supplementary plotting scripts.

Twelve scripts completed successfully. Rewritten CSV outputs matched archived values at the recorded comparison tolerance. Copied-but-unmodified CSVs do not count as recalculation evidence. Logs and machine-readable statuses are written in the new output directory; existing output directories are never overwritten. Image byte identity is not asserted because fonts/rendering can vary.

**Figure 3 limitation:** plotting loads archived Jaccard values; it does not regenerate them under the current missing-reference rule. Inspect both policies with:

```bash
python tools/diagnose_reference_filter.py --input data/processed/merged_results.tsv --output reproduced/reference-policy
```

## Selected plots only

```bash
python tools/reproduce_figures.py --figure 3
python tools/reproduce_figures.py --figure 4
```

These commands use included derived tables. Historical per-figure README commands retain old paths and are not the portable execution interface.

## Scope

- The main TSV has 15,433 rows: 15,145 labeled main-cohort rows and 288 unlabeled rows. Use explicit cohort filters and the separate 288-row FOXL2 TSV; do not blindly concatenate them.
- CancerSIGVAR/hearing-loss and VIP-HL/cancer have no predictions, so no raw file was assigned to those combinations.
- Reference classes/criteria and coordinate mappings are included. Original tool commands, database versions and source-workbook history remain incompletely documented.
- Raw content agreement does not prove historical run identity when alternatives are equivalent.
- Archived values were not changed to resolve the documented ClinGen Jaccard discrepancy.
