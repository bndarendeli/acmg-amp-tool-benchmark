# ACMG Benchmark Analysis Scripts

This repository contains the analysis scripts and outputs for the ACMG/AMP variant interpretation benchmark study.

## Overview

This benchmark evaluates 14 automated ACMG/AMP variant interpretation tools across multiple clinical variant datasets, analyzing both classification-level and evidence-level concordance with expert-curated reference standards.

## Repository Structure

```
.
├── fig2/               # Figure 2: Tool Call Rates & Classification Concordance
├── fig3/               # Figure 3: Criterion-Level Jaccard Analysis
├── fig4/               # Figure 4: Classification Performance Metrics
├── fig5/               # Figure 5: Evidence-Level Concordance Analysis
├── figures/            # Final manuscript figures
├── requirements.txt    # Python dependencies
└── README.md          # This file
```

## Figures

### Figure 2: Tool Call Rates and Classification Concordance

- **Scripts:** `fig2/create_figure2_comprehensive.py`
- **Outputs:** Classification availability and unresolved output patterns
- **Key Finding:** Tool call rates varied substantially across cohorts, with distinct patterns of explicit unresolved and non-returned classifications.

### Figure 3: Criterion-Level Jaccard Analysis

- **Scripts:** `fig3/generate_clingen_jaccard_heatmaps.py`, `fig3/create_figure3_criterion_counts_and_jaccard.py`
- **Outputs:** Jaccard similarity matrices for 28 ACMG/AMP criteria
- **Key Finding:** Evidence-evaluable population: **11,152 ClinGen variants** (257 variants with unavailable reference criteria excluded)
- **Note:** Uses corrected evidence semantics where `'.'` represents missing/unavailable data, not empty evidence sets

### Figure 4: Classification Performance

- **Scripts:** `fig4/create_figure4_dotplot.py`
- **Outputs:** Support-weighted F1 and class-specific recall for P/LP, VUS and B/LB across four cohorts.
- **Key Finding:** Performance varies by dataset and classification type

### Figure 5: Evidence-Level Concordance

- **Scripts:** `fig5/build_fig5_master_table.py`, `fig5/generate_figure5_refined.py`
- **Outputs:** Capability-adjusted evidence concordance, criterion-level decomposition
- **Key Finding:** Raw vs. capability-adjusted Jaccard similarity reveals operational scope limitations

## Installation

```bash
# Clone repository
git clone <repository-url>
cd acmg-benchmark-analysis

# Install dependencies
pip install -r requirements.txt
```

## Requirements

- Python 3.8 or higher
- See `requirements.txt` for package dependencies

## Data Requirements

The analysis scripts expect the following input data structure:

```
data/
├── ground_truth/
│   ├── clingen_28012026_ground_truth.xlsx
│   └── foxl2_ground_truth.xlsx
└── tool_results/
    ├── InterVar_2018/
    ├── BIAS/
    ├── Genebe/
    └── ... (other tools)
```

**Note:** Tool result files and HGMD-derived datasets are not included in this public release pending redistribution rights review. Contact the authors for data access.

## Usage

### Generate Figure 2

```bash
cd fig2
python create_figure2_comprehensive.py
```

### Generate Figure 3

```bash
cd fig3
python generate_clingen_jaccard_heatmaps.py
python create_figure3_criterion_counts_and_jaccard.py
```

### Generate Figure 4

```bash
cd fig4
python create_figure4_dotplot.py
```

### Generate Figure 5

```bash
cd fig5
python build_fig5_master_table.py
python perform_capability_adjusted_analysis.py
python generate_figure5_refined.py
```

## Key Concepts

### Evidence Semantics

**Critical:** Reference criterion field `'.'` represents **missing/unavailable** data, NOT an empty evidence set.

- **Evidence-evaluable variants:** Variants with available reference criteria (not `'.'`)
- **ClinGen evidence-evaluable population:** 11,152 variants (257 VUS variants with `'.'` excluded)
- **FOXL2 evidence-evaluable population:** 288 variants (all have available criteria)

### Capability-Adjusted Analysis

Tools are evaluated only on criteria they can operationally assess given benchmark input constraints:

- **Outside operational scope:** Criterion requires unavailable input (e.g., phenotype, family data)
- **Operationally evaluable but not assigned:** Tool can assess criterion but did not assign it

This distinction separates operational limitations from genuine disagreements.

### 28 Canonical ACMG/AMP Criteria

Analysis focuses on the 28 canonical criteria from the 2015 ACMG/AMP guidelines:

- Very Strong: PVS1
- Strong: PS1-PS4
- Moderate: PM1-PM6
- Supporting: PP1-PP5
- Stand-alone: BA1
- Strong: BS1-BS4
- Supporting: BP1-BP7

Strength modifiers (e.g., `PM2_Supporting`) are collapsed to canonical forms (e.g., `PM2`).

## Outputs

### CSV Data Files

- Jaccard similarity matrices
- Capability-adjusted concordance metrics
- Criterion-level decomposition tables
- Tool performance metrics

### Figure Files

- PNG (300 DPI) for presentations
- PDF (vector) for publication
- SVG (vector) for editing

## Citation

If you use these scripts or data, please cite:

[Citation to be added upon publication]

## License

[License to be determined]

## Contact

For questions about the analysis or data access:

- [Contact information to be added]

## Acknowledgments

This work was supported by [funding information to be added].

## Version History

- **v1.0.0** (2026-09-12): Initial public release
  - Corrected Figure 3 evidence semantics (11,152 ClinGen variants)
  - Renamed Figure 4 → Figure 5 (evidence-level concordance)
  - Excluded archival and development artifacts
  - Clean public-facing documentation

## Notes

- **Archival repository:** The complete development history, including audit trails and reproduction artifacts, is maintained separately
- **HGMD data:** Tool results and datasets containing HGMD-derived information are not included pending redistribution rights review
- **Figure numbering:** This release uses the final manuscript figure numbering (Figures 2-5)
