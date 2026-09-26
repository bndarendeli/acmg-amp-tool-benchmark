# Tool results in the manuscript package

Two frozen processed TSVs and 52 matching raw-output candidates are included. Raw files are stored as gzip streams under `data/tool_results/<tool>/<cohort>/`; original and stored checksums are in `metadata/input_manifest.json`.

Selection followed comparison of 11 merged-input candidates and 153 raw result files, including FOXL2 `results` and `results_fix`. Backups and unrelated cohorts were not bulk-imported. Fifty-one tool/cohort combinations matched parsed classifications and canonical criteria; FOXL2 VIP-HL additionally matched after treating eight Unknown versus missing values as no-calls. Two all-no-call combinations require no source output.

Reference fields are embedded in the processed TSVs. Tool executable versions, original run commands and database snapshots are not reconstructed by packaging. Populate that provenance from actual run records before claiming reproduction of the original tool runs.

See `INPUT_VERIFICATION.md` for the ClinGen Jaccard rule discrepancy requiring reconciliation before publication. No remote repository was created and no files were pushed.
