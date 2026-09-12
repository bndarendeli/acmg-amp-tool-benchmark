# FOXL2 Supplementary Analysis Report

## Purpose

Extend the evidence-concordance analysis from Main Figure 4 (ClinGen) to the
independent FOXL2 dataset to determine whether the observed relationship between
classification concordance and evidence concordance generalizes.

## Pipeline Consistency

✓ Same criterion parser (utils_criteria_parser.py)
✓ Same canonical 28-criterion vocabulary
✓ Same missing-reference handling ('.' = missing/unavailable)
✓ Same capability matrix (tool_acmg_capabilities.csv)
✓ Same Jaccard calculation logic
✓ Same capability-adjusted analysis

## FOXL2 Cohort Flow

**Total FOXL2 variants:** 288
- Variants with documented reference criteria: 288
- Variants without documented reference criteria: 0

**Total variant-tool pairs:** 2191
**Classification-concordant pairs:** 388
**Evidence-evaluable pairs:** 388
  (classification-concordant with documented reference criteria)

## Tool-Level Summary

| Tool | Class Conc | Raw Med J | Adj Med J | ΔJ | Raw Exact | Adj Exact | Op Crit | n Raw | n Adj |
|------|-----------|-----------|-----------|-----|-----------|-----------|---------|-------|-------|
| AutoGVP | 0.269 | 0.500 | 0.750 | +0.250 | 0.000 | 0.143 | 16/28 | 35 | 35 |
| Franklin | 0.347 | 0.500 | 0.500 | +0.000 | 0.000 | 0.000 | 19/28 | 100 | 100 |
| InterVar_2018 | 0.082 | 0.500 | 0.600 | +0.100 | 0.000 | 0.125 | 18/28 | 8 | 8 |
| TAPES | 0.035 | 0.500 | 0.600 | +0.100 | 0.000 | 0.100 | 20/28 | 10 | 10 |
| VIP-HL | 0.021 | 0.417 | 0.875 | +0.458 | 0.000 | 0.500 | 13/28 | 6 | 6 |
| BIAS | 0.205 | 0.333 | 0.500 | +0.167 | 0.000 | 0.186 | 19/28 | 59 | 59 |
| Genebe | 0.413 | 0.333 | 0.600 | +0.267 | 0.000 | 0.084 | 17/28 | 119 | 119 |
| DiabloACMG | 0.279 | 0.286 | 0.400 | +0.114 | 0.000 | 0.103 | 19/28 | 39 | 39 |
| InterVar_2025 | 0.092 | 0.111 | 0.143 | +0.032 | 0.000 | 0.000 | 18/28 | 12 | 12 |

## Overall FOXL2 Statistics

**Mean change in median Jaccard after capability adjustment:** +0.165
**Mean change in exact-match rate after capability adjustment:** +0.138

## Criterion-Level Decomposition

Criteria with missing reference assignments decomposed into:
- Outside operational scope (tool cannot evaluate under benchmark inputs)
- Supported but not assigned (tool can evaluate but did not assign)

| Criterion | Total Missing | Outside Scope | Supported Not Assigned |
|-----------|---------------|---------------|------------------------|

**Scope-dominated (≥80% outside scope):**
  None

**Assignment-dominated (≤20% outside scope):**
  None

**Mixed:**
  None

## Sample Size Assessment

**Total evidence-evaluable predictions across all tools:** 388

**Tools with n ≥ 10:** 7
  - AutoGVP: n=35
  - Franklin: n=100
  - TAPES: n=10
  - BIAS: n=59
  - Genebe: n=119
  - DiabloACMG: n=39
  - InterVar_2025: n=12

**Tools with n < 10 (interpret with caution):** 2
  - InterVar_2018: n=8
  - VIP-HL: n=6

## Comparison to ClinGen Pattern

**ClinGen Main Finding:**
- High classification concordance does not guarantee high evidence concordance
- Capability adjustment reveals operational scope differences
- Missing criteria decompose into scope vs. assignment mechanisms

**FOXL2 Assessment:** Pattern DIFFERS from ClinGen or insufficient variation

## Data Verification

✓ FOXL2 uses same criterion parser as ClinGen
✓ Missing reference criteria ('.' etc.) correctly excluded
✓ Same 28 canonical criteria vocabulary
✓ Same capability matrix and operational definitions
✓ Classification concordance uses exact five-class matching
✓ Evidence-evaluable = classification-concordant + documented reference criteria

## Output Files

**Master Table:**
- `foxl2_fig4_master_table.csv`

**Tool Summary:**
- `foxl2_tool_evidence_summary.csv`

**Criterion Decomposition:**
- `foxl2_criterion_missing_decomposition.csv`

**Audit:**
- `foxl2_missing_reference_audit.txt`

---

*Analysis completed: 2026-09-10*