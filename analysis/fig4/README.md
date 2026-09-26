# Figure 4: Classification Performance (Dot Plot Design)

## Overview

Figure 4 presents classification performance metrics across four cohorts using horizontal multi-marker dot plots with corrected 3-class performance calculations.

## Directory Structure

```
fig4/
├── create_figure4_dotplot.py              # Main figure generation script
├── scripts/
│   └── recalculate_correct_metrics.py    # Metric calculation from confusion matrices
├── data/
│   └── corrected/                         # Corrected metric CSV files
│       ├── ClinGen_corrected_metrics.csv
│       ├── FOXL2_corrected_metrics.csv
│       ├── HGMD+ClinVar_HL_corrected_metrics.csv
│       └── HGMD+ClinVar_Cancer_corrected_metrics.csv
├── figures/                               # Output figures
│   ├── Figure4_Performance_DotPlot.png   (300 DPI)
│   └── Figure4_Performance_DotPlot.pdf   (vector)
├── FIGURE4_DOTPLOT_DESIGN.md             # Design documentation
└── README.md
```

## Quick Start

### Generate Figure 4

```bash
cd /mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_scripts/fig4
python3 create_figure4_dotplot.py
```

### Recalculate Metrics (if needed)

```bash
cd scripts
python3 recalculate_correct_metrics.py
cd ..
python3 create_figure4_dotplot.py
```

## Figure Design

### Visual Hierarchy

**Primary Metric:**
- **Weighted F1 (3-class):** Deep indigo circle, size 75 - visually dominant

**Secondary Metrics:**
- **P/LP Recall:** Vermillion/red-orange triangle, size 38
- **B/LB Recall:** Steel blue square, size 38
- **VUS Recall:** Amber/gold diamond, size 38

### Cohorts (2×2 Layout)

- **Panel a:** ClinGen Master Cohort (n=11,409)
- **Panel b:** FOXL2 Disease-Specific Cohort (n=288)
- **Panel c:** HGMD–ClinVar Hearing Loss (n=1,948)
- **Panel d:** HGMD–ClinVar Cancer Predisposition (n=1,788)

## Corrected Metrics

### Key Corrections

1. **Weighted F1 (3-class):**
   - Calculated from true ground-truth class distribution
   - No-calls count as false negatives for their true class
   - Uses actual full-cohort class support as weights

2. **Recall Metrics:**
   - P/LP Recall: Pathogenic + Likely Pathogenic
   - B/LB Recall: Benign + Likely Benign
   - VUS Recall: VUS (Uncertain Significance)

3. **Call Rate:**
   - Removed from figure (covered in Figure 2)
   - Only valid 5-tier ACMG classifications counted

## Performance Summary

| Cohort | Tools | Mean F1 | Top Performer |
|--------|-------|---------|---------------|
| ClinGen | 14 | 0.587 | AutoGVP (0.887) |
| FOXL2 | 12 | 0.339 | Genebe (0.932) |
| HGMD–ClinVar HL | 13 | 0.477 | TAPES (0.723) |
| HGMD–ClinVar Cancer | 13 | 0.502 | CancerSIGVAR (0.942) |

## Design Features

✅ Horizontal multi-marker dot plots  
✅ Clear visual hierarchy (F1 dominant)  
✅ Colorblind-friendly palette  
✅ Distinct marker shapes  
✅ Tools sorted by Weighted F1 (highest at top)  
✅ Subtle x-axis grid only  
✅ White background  
✅ Nature Communications typography  
✅ 300 DPI publication quality  

## Color Palette

| Metric | Color | Code |
|--------|-------|------|
| Weighted F1 (3-class) | Deep Indigo/Purple | #4B0082 |
| P/LP Recall | Vermillion/Red-Orange | #D94F30 |
| B/LB Recall | Steel Blue | #4682B4 |
| VUS Recall | Amber/Gold | #DAA520 |

## Output Files

- `figures/Figure4_Performance_DotPlot.png` - 300 DPI raster
- `figures/Figure4_Performance_DotPlot.pdf` - Vector (scalable)

## Dependencies

- Python 3.x
- pandas
- numpy
- matplotlib

## Notes

- All metrics calculated from confusion matrices
- Artifact rows (analysis, Genebe_Gene) removed
- Only actual benchmark tool instances displayed
- Panel titles use "HGMD–ClinVar" (not "HGMD+ClinVar")

---

*Last Updated: 2026-09-10*  
*Status: Publication-ready*
