# Input verification and manuscript consistency

## Selected processed inputs

The selected files are unchanged copies of `analysis_last/data/merged_results.tsv` and `analysis_last/data/foxl2_merged_results.tsv`, now in `data/processed/`.

| Cohort | Evaluated rows | Figure 4 metric cells matching at archived four-decimal precision |
| --- | ---: | ---: |
| ClinGen | 11,409 | 154 / 154 |
| HGMD–ClinVar hearing loss | 1,948 | 154 / 154 |
| HGMD–ClinVar cancer | 1,788 | 154 / 154 |
| FOXL2 | 288 | 132 / 132 |

All 594 metric cells matched, including undefined values in zero-support classes. All 54 Figure 2 call-count records matched using both the archived non-null/non-Unknown rule and strict valid-class normalization.

ClinGen and FOXL2 Figure 5 master-table keys, reference/tool classifications and canonical criteria matched without extra or missing keys. An isolated run regenerated the full master table, capability-adjusted analyses and supplementary outputs; rewritten CSVs matched the archive. All 12 executed analysis/plotting scripts returned exit status zero. Unchanged copied CSVs are not evidence of recalculation.

## ClinGen Figure 3 discrepancy requiring author review

The archived criterion-level Jaccard table is reproduced exactly (392/392 intersection/union pairs) when all 11,409 ClinGen rows are included and missing reference criteria are treated as empty sets. Applying the current script's exclusion of 257 missing-reference rows matches only 283/392 pairs: 109 count pairs and 94 Jaccard values differ (absolute tolerance 1e-12 for coefficient comparisons).

The current script and `EVIDENCE_SEMANTICS_CONSISTENCY_REPORT.md` describe exclusion of missing reference criteria. The archived table therefore does not implement the documented rule. Auditing the historical reference-workbook join gave the same discrepancy, with 11,409 unique variants and no join expansion.

FOXL2 criterion-level Jaccard matched 336/336 pairs with the current rule. Figure 3 plotting completed using archived tables; successful plotting does not resolve the ClinGen discrepancy.

No archived table or manuscript value was overwritten. `metadata/audit/reference-policy-01/policy_comparison.csv` contains both calculations for review. The author should choose the scientifically intended definition, then align the calculation, manuscript text, affected figure values and method notes together.

## Why alternatives were not selected

- Older merged files contain `Uncertain_significance` in hearing-loss/cancer reference labels, which the current strict normalization function does not accept. Those candidates cannot reproduce the current pipeline unchanged.
- `merged_results_fixed.tsv` agrees for the three labeled cohorts but differs from the selected file in 288 reference-label cells on rows without a Dataset label.
- The main merged snapshot has 288 unlabeled rows. The analyses filter the three main cohorts and use the separate FOXL2 file, avoiding double counting.
- FOXL2 outputs in `results_fix` support the selected merged table; earlier `results` versions are not interchangeable.

## Raw-output traceability

Of 153 examined files, candidates were selected for 52 populated tool/cohort combinations. Fifty-one matched parsed classifications and canonical criterion fields directly. FOXL2 VIP-HL had eight parsed Unknown versus missing merged classifications and matching criteria. These are equivalent no-calls for the verified metrics, but not identical encodings. Its raw file has two duplicate variant keys, handled by keeping the first as in the historical merge code.

CancerSIGVAR/hearing-loss and VIP-HL/cancer contain only no-calls and empty criteria. No irrelevant file was selected merely to fill those combinations. Equivalent candidates are recorded in metadata; content matching is not proof of historical run provenance.

## Evidence

- `metadata/input_manifest.json`: selected paths, original/stored hashes, sizes and selection status.
- `metadata/selected_tool_results.csv`: raw-file/cohort mapping.
- `metadata/audit/input-audit-01/`: 11-input inventory and quantitative comparisons.
- `metadata/audit/raw-audit-01/`: full raw candidate inventory and comparisons.
- `metadata/audit/selected-details-01/`: alternative-input differences and VIP-HL records.
- `metadata/audit/full-reproduction-01/`: execution statuses and CSV comparisons.

Original source directories remain unchanged. Execution products are isolated under ignored `reproduced/`.
