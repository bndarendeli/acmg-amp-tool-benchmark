# Figure and supplementary-table guide

Mappings use the supplied draft's subject matter and existing filenames; they do not certify pixel-identical manuscript versions.

| Manuscript item / analysis | Tables | Available figure exports |
| --- | --- | --- |
| Figure 2: classification availability | `analysis/fig2/corrected_call_rates.csv` | `analysis/fig2/figures/Figure2_Comprehensive.{pdf,png}` |
| Figure 3: ClinGen criterion counts and Jaccard | `analysis/fig3/figures/figure3_clingen_criterion_counts.csv`, `figure3_clingen_jaccard.csv`; `analysis/fig3/jaccard_per_criterion_clingen_28012026.csv` | `analysis/fig3/figures/Figure3_ClinGen_Counts_Jaccard.{pdf,png}` |
| Figure 4: classification performance | Four CSVs in `analysis/fig4/data/corrected/` | `analysis/fig4/figures/Figure4_Performance_DotPlot.{pdf,png}` |
| Figure 5: evidence concordance | `analysis/fig5/outputs/fig4_master_variant_tool.csv`, `fig4_evidence_distributions.csv`, `fig4_capability_adjusted_summary.csv`, `fig4_missing_criteria_decomposition.csv` | Multiple historical variants in `analysis/fig5/outputs/figures/`, including `Figure5_redesign_preview_v2.pdf` and `fig4_final.pdf`; final version selection remains with the author |
| Supplementary S1: InterVar comparison | Requires external merged predictions to rebuild | `analysis/fig2/figures/Supp_Fig_InterVar_Comparison.{pdf,png}` |
| Supplementary S2: FOXL2 criterion counts and Jaccard | `analysis/fig3/figures/figure3_foxl2_criterion_counts.csv`, `figure3_foxl2_jaccard.csv` | `analysis/fig3/figures/Supp_Fig2_FOXL2_Counts_Jaccard.pdf` |
| Supplementary S3: FOXL2 evidence concordance | `analysis/fig5/outputs/foxl2_supplementary/` | `analysis/fig5/outputs/figures/Supp_Fig3_FOXL2_Supplementary_Evidence_Concordance.pdf` |
| Supplementary tables S3–S5: tool scope and ClinGen evidence | Three CSVs in `analysis/fig5/outputs/clingen_supplementary/` | Not applicable |
| Supplementary tables S6–S7: FOXL2 evidence | `analysis/fig5/outputs/foxl2_supplementary/foxl2_capability_adjusted_summary.csv`, `foxl2_criterion_missing_decomposition.csv` | Not applicable |

No Figure 1 source file was present. Older exports such as `Figure4_Performance_2x2.pdf`, individual SVG panels and preview files remain in the archive for provenance. Filenames containing “final” are historical labels, not a packaging-time selection.
