# Validation report — 2026-10-28

The repository was validated and its complete selected-figure pipeline was run in the Ubuntu WSL `acmg-benchmarking` Conda environment (Python 3.9). The run used `python tools/run_all.py --from-raw --output /tmp/acmg-benchmark-fig6-chain-check`; it created an isolated output copy and did not overwrite the repository's supplied figures or data.

- Repository validation passed: 32 manifested inputs, 42 Python files parsed, and expected cohort sizes were found.
- All 22 pipeline steps, including raw parsing and the Figure 6 builders, completed successfully.
- The regenerated `merged_results_noHGMD.tsv` was byte-identical to the bundled table (SHA-256 recorded in `validation_results.json`).
- Six regenerated Fig6 input tables matched the supplied tables, including ClinGen and FOXL2 master/summary tables, evidence distributions, total Jaccard tables, and the capability matrix.
- The evidence distribution builder uses classification-concordant, evidence-evaluable rows and reproduces the bundled interval-bin values.
- Figure 6 inputs come from the current ClinGen and FOXL2 files only. The combined no-HGMD TSV's FOXL2 rows are excluded from the ClinGen input and the dedicated FOXL2 file is used to avoid duplicate records.
- Original annotation tools and external online services are not rerun. Plot appearance/PDF bytes can vary across plotting-library and font versions.

Run `python tools/validate_repository.py` to verify bundled input and provenance hashes and portable script syntax. Run `python tools/run_all.py --from-raw` to reproduce all figures and tables in a new output copy.
