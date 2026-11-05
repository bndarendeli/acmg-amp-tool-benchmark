# Figure 6 and Supplementary Figure S4

The repository includes the selected Figure 6 and Supplementary Figure S4, their plotting scripts, and the analysis scripts and reference table needed to rebuild their input tables from bundled ClinGen and FOXL2 predictions.

| Figure | Script | Outputs |
| --- | --- | --- |
| Figure 6 (selected five-panel v3 design) | `generate_figure6.py` | `outputs/figures/Figure6_five_panel_preview_v3.pdf` |
| Supplementary Figure S4 (FOXL2) | `generate_supplementary_figure_s4.py` | `outputs/foxl2_supplementary/Supplementary_Figure_S4.pdf` |

Run the complete repository pipeline from the repository root with:

```bash
python tools/run_all.py
```

To rebuild merged predictions from the bundled raw tool outputs first, use `python tools/run_all.py --from-raw`.

The upstream scripts run in this order: build the variant/tool master from `data/processed/merged_results_noHGMD.tsv` and `data/processed/foxl2_merged_results.tsv`; build the capability matrix using the bundled `reference/tool_acmg_capabilities.csv`; calculate capability-adjusted ClinGen and FOXL2 summaries; build evidence distributions; then calculate total criterion Jaccard tables.

Figure 6 reads five tables in `outputs/`: `fig6_evidence_distributions.csv`,
`fig6_master_variant_tool.csv`, `fig6_capability_adjusted_summary.csv`,
`clingen_total_criterion_jaccard.csv`, and `tool_criterion_capability_matrix.csv`.
Supplementary Figure S4 reads `outputs/foxl2_supplementary/foxl2_total_criterion_jaccard.csv`.
These tables are rebuilt by the upstream scripts above. The master table and adjusted summary are also used by Figures 3 and 5. The analysis uses only ClinGen and FOXL2 benchmark data; no HGMD/ClinVar merged dataset is an input.
