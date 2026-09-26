# Data

`processed/` contains frozen tool predictions, reference classifications and reference ACMG criteria. `tool_results/<tool>/<cohort>/` contains losslessly compressed matching raw-output candidates. Original filenames have `.gz` appended.

See `metadata/input_manifest.json` for hashes, provenance and verification status. `python tools/rebuild_tool_predictions.py --output reproduced/raw-rerun` decompresses/reparses files in temporary directories and compares normalized predictions with the snapshots.

The main merged TSV contains 288 rows without a Dataset label. Filter the three labeled main cohorts and use the separate FOXL2 TSV for FOXL2; do not blindly concatenate both files. Reference fields and genome-build mappings are embedded in these tables.

These are benchmark records, not tool executables or complete annotation databases. No new redistribution license is implied.
