# Figure 3: Unified ACMG Criteria Heatmap

This directory contains all scripts and data files needed to generate Figure 3 (Unified ACMG Heatmap - Two Panel Design).

## Output Figure

- **Figure**: `Figure3_Unified_ACMG_Heatmap_TwoPanel.png` / `.pdf`
- **Location**: `figures/` (analysis/fig3/figures/)

## Files

### Scripts

#### Figure Generation

- `create_figure3_unified_acmg_heatmap.py` - Main script to generate the two-panel ACMG heatmap

#### Data Generation (Optional - data files already included)

- `generate_clingen_jaccard_heatmaps.py` - Generates Jaccard similarity data for ClinGen dataset
- `analyze_foxl2_comprehensive.py` - Generates Jaccard similarity data for FOXL2 dataset

### Data Files (Pre-generated)

- `jaccard_per_criterion_clingen_28012026.csv` - Jaccard similarity data for ClinGen dataset (Panel A)
- `jaccard_per_criterion_foxl2.csv` - Jaccard similarity data for FOXL2 dataset (Panel B)

## Figure Description

**Panel A (Left)**: ClinGen Master Cohort

- Tool alignment matrix showing Jaccard similarity coefficients
- 14 tools × ACMG criteria
- Criteria grouped by categories: Population, Computational/Predictive, Functional

**Panel B (Right)**: FOXL2 Dataset

- Same tool alignment matrix for FOXL2 validation cohort
- Same tools and criteria organization as Panel A

## How to Run

### Quick Start (Using Pre-generated Data)

```bash
cd /mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_scripts/fig3
python3 create_figure3_unified_acmg_heatmap.py
```

### Full Pipeline (Regenerate Data from Scratch)

If you need to regenerate the Jaccard similarity data:

1. **Generate ClinGen data:**

```bash
python3 generate_clingen_jaccard_heatmaps.py
```

This will create `jaccard_per_criterion_clingen_28012026.csv`

2. **Generate FOXL2 data:**

```bash
python3 analyze_foxl2_comprehensive.py
```

This will create `jaccard_per_criterion_foxl2.csv` (as `foxl2_jaccard_per_criterion.csv`)

3. **Generate the figure:**

```bash
python3 create_figure3_unified_acmg_heatmap.py
```

## Requirements

- Python 3.x
- matplotlib
- seaborn
- numpy
- pandas

## Output

The script will generate:
- `Figure3_Unified_ACMG_Heatmap_TwoPanel.png` (300 DPI, ~552 KB)
- `Figure3_Unified_ACMG_Heatmap_TwoPanel.pdf` (vector, ~48 KB)

**Output Location**: `figures/` (analysis/fig3/figures/)

## Notes

- The script uses relative paths (`base_dir = Path(__file__).parent.parent`)
- Data files must be in the same directory as the script
- Output directory will be created automatically if it doesn't exist
