# ACMG/AMP tool benchmark

Analysis data, scripts and figures for the ClinGen and FOXL2 benchmark.
This release contains **11,409 ClinGen variants** and **288 FOXL2 variants**.

## Install

Clone the repository into a directory with enough space for a separate analysis copy:

```bash
git clone https://github.com/bndarendeli/acmg-amp-tool-benchmark.git
cd acmg-amp-tool-benchmark
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS/WSL, or
`.venv\Scripts\Activate.ps1` in Windows PowerShell. Then install:

```bash
python -m pip install -r requirements.txt
python tools/validate_repository.py
```

Python 3.9 or later is required. The end-to-end run was performed with the `acmg-benchmarking` Conda environment in Ubuntu WSL.
No personal drive paths, external source directories or Conda environment names are required.

## Reproduce the analyses

```bash
python tools/run_all.py
```

The runner copies the bundled inputs and scripts into a new `reproduced/<timestamp>/`
directory and executes the current analyses in dependency order. Each script has a log
under that copy's `logs/` directory. Supplied manuscript outputs remain unchanged.
Figures 2–6 and the selected supplementary figures are generated under their respective
copied figure directories.

To rebuild merged predictions from the bundled tool outputs first:

```bash
python tools/run_all.py --from-raw
```

To choose another new output directory:

```bash
python tools/run_all.py --output reproduced/my-run
```

Individual scripts can also be called by path from any working directory; they write
into their own repository figure directory. Use the runner to preserve supplied outputs.

## Contents

| Directory | Contents |
| --- | --- |
| `analysis/fig2` | Coverage, criteria, parsers and InterVar comparison |
| `analysis/fig3` | Criterion counts and Jaccard heatmaps |
| `analysis/fig4` | Two-cohort classification metrics and figures |
| `analysis/fig5` | Selected Pearson correlation heatmap and its inputs |
| `analysis/fig6` | Figure 6, Supplementary Figure S4 and their inputs |
| `data/processed` | Current merged predictions and the separate FOXL2 analysis input |
| `data/reference` | ClinGen-only reference rows and the FOXL2 reference workbook |
| `data/tool_results` | 26 compressed original outputs: 14 ClinGen and 12 FOXL2 |
| `data/supporting` | FOXL2 classification summary used by the call-rate script |
| `metadata` | Input checksums, source-copy provenance and verification report |
| `tools` | Portable runner, input unpacking, validation and optional WSL logging |

## Reproduction scope

The pipeline rebuilds the selected manuscript figures and their analysis tables from bundled
ClinGen and FOXL2 inputs. The Figure 6 chain now builds the master variant/tool table, capability
matrix, capability-adjusted summaries, evidence distributions and total criterion Jaccard tables.
It reads ClinGen predictions from `data/processed/merged_results_noHGMD.tsv` and FOXL2 predictions
from `data/processed/foxl2_merged_results.tsv`; the FOXL2 rows in the combined file are filtered
out to prevent duplication.

Each figure folder contains the selected manuscript figure files, their producing scripts, and required data. Follow this README
for installation and execution. There is no shared top-level `figures` directory.
Font availability and plotting-library versions may affect appearance or PDF bytes.
See `metadata/validation_report.md` for the actual numerical checks and tested environment.
