# Figure 3 Update Summary

## Overview

Updated Figure 3 to show criterion assignment counts (Panel a) alongside Jaccard similarity (Panel b) in a side-by-side layout matching the original figure's visual style.

## Changes Made

### Layout & Visual Design
- **Two-panel side-by-side layout** matching original Figure 3 proportions
- **Figure size:** 16×10 inches
- **Panel spacing:** wspace=0.3
- **Width ratios:** 1.1:1 (counts panel slightly wider for Ground truth column)
- **Category separators:** Black horizontal lines between Population/Computational/Functional groups
- **Category labels:** Vertical text with brackets on left side (matching original)
- **Panel labels:** Lowercase 'a' and 'b' in bold (fontsize=20)
- **No panel titles** (as requested)

### Panel a: Criterion Assignment Counts
- **Data:** Number of times each criterion was assigned by each tool
- **Columns:** All 14 tools + "Ground truth" (reference counts)
- **Rows:** 18 selected ACMG/AMP criteria (from original Figure 3)
- **Color scale:** `plasma` (purple to orange/yellow)
- **Annotations:** Count values displayed in each cell
- **Colorbar:** None (cleaner look)

### Panel b: Jaccard Similarity
- **Data:** Criterion-level Jaccard similarity (unchanged from original)
- **Columns:** All 14 tools
- **Rows:** Same 18 criteria (consistent ordering)
- **Color scale:** `RdBu_r` (red-blue, preserved from original)
- **Value range:** 0-1
- **Annotations:** Jaccard values (2 decimal places)
- **Colorbar:** Right side with "Jaccard Similarity Coefficient" label

### Selected Criteria (18 total)

**Population (5):**
- BA1, BS1, BS2, PM2, PS4

**Computational and Predictive (9):**
- BP1, BP3, BP4, BP7, PP3, PM4, PM5, PS1, PVS1

**Functional (4):**
- BS3, PP2, PM1, PS3

## Data Processing

### Criterion Counts Calculation
- Uses same criterion parser as Figure 4 (`utils_criteria_parser.py`)
- Treats missing reference values (`.`) as unavailable (not empty sets)
- Counts only canonical criteria (strength modifiers normalized)
- Ground truth counts calculated from unique variants only

### Tool Coverage
- **ClinGen:** 14 tools, 104,152 variant-tool pairs
- **FOXL2:** 12 tools, 2,191 variant-tool pairs

### Quality Control
✓ Same criterion parser and normalization as all other figures
✓ Missing reference handling consistent with Figure 4 pipeline
✓ Criterion ordering identical between panels a and b
✓ Category boundaries aligned across both panels
✓ No data modifications or calculations changed

## Output Files

### ClinGen Main Figure
- `Figure3_ClinGen_Counts_Jaccard.png` (300 dpi)
- `Figure3_ClinGen_Counts_Jaccard.pdf`
- `figure3_clingen_criterion_counts.csv` (data matrix)
- `figure3_clingen_jaccard.csv` (data matrix)

### FOXL2 Supplementary Figure
- `Figure3_FOXL2_Counts_Jaccard.png` (300 dpi)
- `Figure3_FOXL2_Counts_Jaccard.pdf`
- `figure3_foxl2_criterion_counts.csv` (data matrix)
- `figure3_foxl2_jaccard.csv` (data matrix)

**Location:** `/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis/fig3/figures/`

## Script

**Main script:** `create_figure3_criterion_counts_and_jaccard.py`

**Key functions:**
- `calculate_criterion_counts()` - Count criterion assignments per tool
- `load_jaccard_matrix()` - Load existing Jaccard data
- `create_figure_for_dataset()` - Generate two-panel figure with category separators
- `save_matrices_as_csv()` - Export data matrices

## Consistency Checks

✓ Uses corrected ClinGen pipeline (Ground_Truth_ACMG='.' handled correctly)
✓ Same 18 criteria as original Figure 3
✓ Category groupings preserved (Population/Computational/Functional)
✓ Visual style matches original figure (separators, brackets, labels)
✓ Panel proportions and spacing match original
✓ No changes to Jaccard calculation or data
✓ Tool ordering consistent with other figures

## Notes

- Panel a colorbar removed for cleaner appearance (counts are annotated)
- Panel b colorbar preserved on right side (standard for Jaccard)
- Category separators and labels applied to both panels
- Ground truth column added only to Panel a (counts)
- Both panels use same Y-axis (criterion) ordering
- Lowercase panel labels per manuscript style guide

---

*Generated: 2026-09-10*
*Script: create_figure3_criterion_counts_and_jaccard.py*
