# Validation results

Validated on 2026-09-12 using Ubuntu WSL2, Conda base, Python 3.12.11.

- 125 archived analysis files and 33 archived CSV tables; original scientific values retained.
- 54 added input files: 52 gzip-compressed raw-output candidates and two merged TSVs.
- Added inputs total 47.93 MiB; largest stored input is 11.08 MiB.
- Gzip decompression was checked against original source SHA-256 hashes during compression.
- All 52 populated tool/cohort combinations were reparsed from the packaged files and matched normalized five-tier classifications and canonical criteria. Two all-no-call combinations were separately verified. Eight FOXL2 VIP-HL Unknown/missing encoding differences remain documented.
- All 594 Figure 4 metric cells and 54 Figure 2 call counts matched the selected inputs.
- Twelve analysis/plotting scripts executed successfully in an isolated copy. Rewritten CSV outputs matched the archive; images are not claimed to be byte-identical.
- ClinGen Figure 3: excluding 257 missing-reference variants changes 109 intersection/union count pairs and 94 Jaccard values. The archived table instead matches inclusion of those rows. This inconsistency is preserved for author review.
- No original source files were edited; no remote repository, commit or push was created.

Run `python tools/validate_repository.py` for current integrity checks. See `INPUT_VERIFICATION.md` for comparison scope and limitations, and `metadata/audit/` for machine-readable evidence.
