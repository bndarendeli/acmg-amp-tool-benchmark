# Archive notes

Prepared on 2026-09-12 from the supplied `analysis_scripts` directory. The attached `acmg-article1.pdf` was read for scientific context only and was not copied into the release.

## Provenance

`source_manifest.json` maps each original relative path to its repository path, byte size and SHA-256 hash. Data and figure contents are unchanged; selected script output paths and navigation documents were updated. The manifest records both original and packaged SHA-256 hashes. Compiled Python caches are excluded. `table_schema.json` and `DATA_CATALOG.md` record observed CSV headers and row counts, not inferred biological definitions.

## Historical inconsistencies

- The PDF first page contains placeholder title/authors. Citation metadata must come from the final article.
- Existing README files contain local machine paths, stale filenames, prior metrics and outdated figure numbering. They are preserved as historical documents, not verified execution instructions. Use the new root documentation to navigate the archive.
- Figure 5 scripts and tables often retain `fig4` names. Some scripts import a parser from a nonexistent `fig1` directory; the available parser is in `fig2`.
- The historical Figure 4 README describes confusion-matrix metrics, while a variant-level recalculation script and corrected CSVs are also present. The old confusion-matrix script is explicitly named `OLD_CONFUSION_MATRIX` and should not be treated as the canonical calculation.
- The Figure 5 README says all FOXL2 reference variants are pathogenic, whereas the supplied manuscript reports five reference VUS. No data were changed to reconcile this statement.

## Release scope

The initial archive was extended after input verification with two merged prediction/reference tables and 52 matching raw-output candidates. These files are recorded separately in `metadata/input_manifest.json`; the original analysis manifest is unchanged. Reference labels and criteria are embedded in the merged TSVs. The package supports raw-output reparsing and isolated Figures 2–5 analysis reproduction, but does not rerun the original tools or recover all original database/run versions. See `INPUT_VERIFICATION.md` for the confirmed ClinGen Jaccard policy discrepancy and other selection details.

Licensing and final article metadata remain unspecified. No LICENSE or CITATION.cff with invented author or publication details has been added. Packaging checks do not establish redistribution rights for third-party derived data or scientific agreement with every manuscript value.
