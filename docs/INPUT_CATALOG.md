# Processed and raw input catalog

Two processed TSVs contain reference data and predictions. Counts below describe physical records, not automatically the evaluated population.

## [merged_results.tsv](../data/processed/merged_results.tsv)

Rows: 15,433; columns: 40.

Columns: `Variant_Key`, `Variant_Key_hg19`, `Chr`, `Pos`, `Chr_hg19`, `Pos_hg19`, `Ref`, `Alt`, `Ground_Truth_Classification`, `Ground_Truth_ACMG`, `Dataset`, `InterVar_2018_Classification`, `InterVar_2018_ACMG_Criteria`, `InterVar_2025_Classification`, `InterVar_2025_ACMG_Criteria`, `BIAS_Classification`, `BIAS_ACMG_Criteria`, `CharGer_Local_Classification`, `CharGer_Local_ACMG_Criteria`, `CharGer_Online_Classification`, `CharGer_Online_ACMG_Criteria`, `DiabloACMG_Classification`, `DiabloACMG_ACMG_Criteria`, `Exomiser_Classification`, `Exomiser_ACMG_Criteria`, `Genebe_Classification`, `Genebe_ACMG_Criteria`, `TAPES_Classification`, `TAPES_ACMG_Criteria`, `Franklin_Classification`, `Franklin_ACMG_Criteria`, `AutoGVP_Classification`, `AutoGVP_ACMG_Criteria`, `VIP-HL_Classification`, `VIP-HL_ACMG_Criteria`, `CancerSIGVAR_Classification`, `CancerSIGVAR_ACMG_Criteria`, `CPSR_Classification`, `CPSR_ACMG_Criteria`, `HGVS`

## [foxl2_merged_results.tsv](../data/processed/foxl2_merged_results.tsv)

Rows: 288; columns: 29.

Columns: `Variant_Key`, `Variant_Key_hg19`, `HGVS`, `Ground_Truth_Classification`, `Ground_Truth_ACMG`, `InterVar_2018_Classification`, `InterVar_2018_ACMG_Criteria`, `InterVar_2025_Classification`, `InterVar_2025_ACMG_Criteria`, `BIAS_Classification`, `BIAS_ACMG_Criteria`, `CharGer_Local_Classification`, `CharGer_Local_ACMG_Criteria`, `CharGer_Online_Classification`, `CharGer_Online_ACMG_Criteria`, `DiabloACMG_Classification`, `DiabloACMG_ACMG_Criteria`, `Exomiser_Classification`, `Exomiser_ACMG_Criteria`, `Genebe_Classification`, `Genebe_ACMG_Criteria`, `TAPES_Classification`, `TAPES_ACMG_Criteria`, `Franklin_Classification`, `Franklin_ACMG_Criteria`, `AutoGVP_Classification`, `AutoGVP_ACMG_Criteria`, `VIP-HL_Classification`, `VIP-HL_ACMG_Criteria`

## Raw tool outputs

| Tool | Cohort | File | Verification |
| --- | --- | --- | --- |
| AutoGVP | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38_FINAL_CLEAN-autogvp-annotated-full_acmg_criteria_parsed.tsv.gz](../data/tool_results/AutoGVP/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38_FINAL_CLEAN-autogvp-annotated-full_acmg_criteria_parsed.tsv.gz) | canonical_prediction_and_criteria_match |
| AutoGVP | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38_FINAL_CLEAN-autogvp-annotated-full_acmg_criteria_parsed.tsv.gz](../data/tool_results/AutoGVP/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38_FINAL_CLEAN-autogvp-annotated-full_acmg_criteria_parsed.tsv.gz) | canonical_prediction_and_criteria_match |
| AutoGVP | clingen_28012026 | [clingen_28012026_hg38_FINAL_CLEAN-autogvp-annotated-full_acmg_criteria_parsed.tsv.gz](../data/tool_results/AutoGVP/clingen_28012026/clingen_28012026_hg38_FINAL_CLEAN-autogvp-annotated-full_acmg_criteria_parsed.tsv.gz) | canonical_prediction_and_criteria_match |
| AutoGVP | foxl2 | [foxl2_hg38_autogvp.tsv.gz](../data/tool_results/AutoGVP/foxl2/foxl2_hg38_autogvp.tsv.gz) | canonical_prediction_and_criteria_match |
| BIAS | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38.tsv.gz](../data/tool_results/BIAS/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38.tsv.gz) | canonical_prediction_and_criteria_match |
| BIAS | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38.tsv.gz](../data/tool_results/BIAS/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38.tsv.gz) | canonical_prediction_and_criteria_match |
| BIAS | clingen_28012026 | [clingen_28012026_hg38.tsv.gz](../data/tool_results/BIAS/clingen_28012026/clingen_28012026_hg38.tsv.gz) | canonical_prediction_and_criteria_match |
| BIAS | foxl2 | [foxl2_hg38_BIAS.tsv.gz](../data/tool_results/BIAS/foxl2/foxl2_hg38_BIAS.tsv.gz) | canonical_prediction_and_criteria_match |
| CPSR | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38_grch38.cpsr.grch38.classification.tsv.gz](../data/tool_results/CPSR/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38_grch38.cpsr.grch38.classification.tsv.gz) | canonical_prediction_and_criteria_match |
| CPSR | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38_grch38.cpsr.grch38.classification.tsv.gz](../data/tool_results/CPSR/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38_grch38.cpsr.grch38.classification.tsv.gz) | canonical_prediction_and_criteria_match |
| CPSR | clingen_28012026 | [clingen_28012026_hg38_grch38.cpsr.grch38.classification.tsv.gz](../data/tool_results/CPSR/clingen_28012026/clingen_28012026_hg38_grch38.cpsr.grch38.classification.tsv.gz) | canonical_prediction_and_criteria_match |
| CancerSIGVAR | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg19_CancerSIGVAR_input_formatted_Cancer_SIGVAR_Results.tsv.gz](../data/tool_results/CancerSIGVAR/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg19_CancerSIGVAR_input_formatted_Cancer_SIGVAR_Results.tsv.gz) | canonical_prediction_and_criteria_match |
| CancerSIGVAR | clingen_28012026 | [clingen_cancerpredisposition_28012026_hg19_CancerSIGVAR_input_formatted_Cancer_SIGVAR_Results.tsv.gz](../data/tool_results/CancerSIGVAR/clingen_28012026/clingen_cancerpredisposition_28012026_hg19_CancerSIGVAR_input_formatted_Cancer_SIGVAR_Results.tsv.gz) | canonical_prediction_and_criteria_match |
| CharGer_Local | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38_CharGer_local_annotated_VEPv97_results.tsv.gz](../data/tool_results/CharGer_Local/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38_CharGer_local_annotated_VEPv97_results.tsv.gz) | canonical_prediction_and_criteria_match |
| CharGer_Local | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38_CharGer_local_annotated_VEPv97_results.tsv.gz](../data/tool_results/CharGer_Local/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38_CharGer_local_annotated_VEPv97_results.tsv.gz) | canonical_prediction_and_criteria_match |
| CharGer_Local | clingen_28012026 | [clingen_28012026_hg38_fixed_CharGer_local_annotated_VEPv97_results.tsv.gz](../data/tool_results/CharGer_Local/clingen_28012026/clingen_28012026_hg38_fixed_CharGer_local_annotated_VEPv97_results.tsv.gz) | canonical_prediction_and_criteria_match |
| CharGer_Local | foxl2 | [foxl2_hg38_CharGer_local_annotated_VEPv97.tsv.gz](../data/tool_results/CharGer_Local/foxl2/foxl2_hg38_CharGer_local_annotated_VEPv97.tsv.gz) | canonical_prediction_and_criteria_match |
| CharGer_Online | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38_CharGer_online_annotated_and_local_ClinVar_results.tsv.gz](../data/tool_results/CharGer_Online/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38_CharGer_online_annotated_and_local_ClinVar_results.tsv.gz) | canonical_prediction_and_criteria_match |
| CharGer_Online | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38_CharGer_online_annotated_and_local_ClinVar_results.tsv.gz](../data/tool_results/CharGer_Online/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38_CharGer_online_annotated_and_local_ClinVar_results.tsv.gz) | canonical_prediction_and_criteria_match |
| CharGer_Online | clingen_28012026 | [clingen_28012026_hg38_fixed_CharGer_online_annotated_and_local_ClinVar_results.tsv.gz](../data/tool_results/CharGer_Online/clingen_28012026/clingen_28012026_hg38_fixed_CharGer_online_annotated_and_local_ClinVar_results.tsv.gz) | canonical_prediction_and_criteria_match |
| CharGer_Online | foxl2 | [foxl2_hg38_CharGer_online.tsv.gz](../data/tool_results/CharGer_Online/foxl2/foxl2_hg38_CharGer_online.tsv.gz) | canonical_prediction_and_criteria_match |
| DiabloACMG | HGMD_Clinvar_Cancer | [HGMD_hg38_Clinvar_Cancer_annotated.tsv.gz](../data/tool_results/DiabloACMG/HGMD_Clinvar_Cancer/HGMD_hg38_Clinvar_Cancer_annotated.tsv.gz) | canonical_prediction_and_criteria_match |
| DiabloACMG | HGMD_Clinvar_HL | [HGMD_hg38_Clinvar_HL_annotated.tsv.gz](../data/tool_results/DiabloACMG/HGMD_Clinvar_HL/HGMD_hg38_Clinvar_HL_annotated.tsv.gz) | canonical_prediction_and_criteria_match |
| DiabloACMG | clingen_28012026 | [clingen_hg38_28012026_annotated.tsv.gz](../data/tool_results/DiabloACMG/clingen_28012026/clingen_hg38_28012026_annotated.tsv.gz) | canonical_prediction_and_criteria_match |
| DiabloACMG | foxl2 | [foxl2_hg38_DiabloACMG.tsv.gz](../data/tool_results/DiabloACMG/foxl2/foxl2_hg38_DiabloACMG.tsv.gz) | canonical_prediction_and_criteria_match |
| Exomiser | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38-exomiser.variants.tsv.gz](../data/tool_results/Exomiser/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38-exomiser.variants.tsv.gz) | canonical_prediction_and_criteria_match |
| Exomiser | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38-exomiser.variants.tsv.gz](../data/tool_results/Exomiser/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38-exomiser.variants.tsv.gz) | canonical_prediction_and_criteria_match |
| Exomiser | clingen_28012026 | [clingen_28012026_hg38-exomiser.variants.tsv.gz](../data/tool_results/Exomiser/clingen_28012026/clingen_28012026_hg38-exomiser.variants.tsv.gz) | canonical_prediction_and_criteria_match |
| Exomiser | foxl2 | [foxl2_hg38_exomiser.tsv.gz](../data/tool_results/Exomiser/foxl2/foxl2_hg38_exomiser.tsv.gz) | canonical_prediction_and_criteria_match |
| Franklin | HGMD_Clinvar_Cancer | [merged_HGMD_Clinvar_Cancer_hg19_last_single_snp_variants.csv.gz](../data/tool_results/Franklin/HGMD_Clinvar_Cancer/merged_HGMD_Clinvar_Cancer_hg19_last_single_snp_variants.csv.gz) | canonical_prediction_and_criteria_match |
| Franklin | HGMD_Clinvar_HL | [merged_HGMD_Clinvar_HL_hg19_last_single_snp_variants.csv.gz](../data/tool_results/Franklin/HGMD_Clinvar_HL/merged_HGMD_Clinvar_HL_hg19_last_single_snp_variants.csv.gz) | canonical_prediction_and_criteria_match |
| Franklin | clingen_28012026 | [merged_all_clingen_28012026_hg19_last_single_snp_variants.csv.gz](../data/tool_results/Franklin/clingen_28012026/merged_all_clingen_28012026_hg19_last_single_snp_variants.csv.gz) | canonical_prediction_and_criteria_match |
| Franklin | foxl2 | [foxl2_hg19_franklin.csv.gz](../data/tool_results/Franklin/foxl2/foxl2_hg19_franklin.csv.gz) | canonical_prediction_and_criteria_match |
| Genebe | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38_genebe.vcf.gz](../data/tool_results/Genebe/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38_genebe.vcf.gz) | canonical_prediction_and_criteria_match |
| Genebe | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38_genebe.vcf.gz](../data/tool_results/Genebe/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38_genebe.vcf.gz) | canonical_prediction_and_criteria_match |
| Genebe | clingen_28012026 | [clingen_28012026_hg38_genebe.vcf.gz](../data/tool_results/Genebe/clingen_28012026/clingen_28012026_hg38_genebe.vcf.gz) | canonical_prediction_and_criteria_match |
| Genebe | foxl2 | [foxl2_hg38_genebe.vcf.gz](../data/tool_results/Genebe/foxl2/foxl2_hg38_genebe.vcf.gz) | canonical_prediction_and_criteria_match |
| InterVar_2018 | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38_fixed.hg38_multianno.txt.gz](../data/tool_results/InterVar_2018/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38_fixed.hg38_multianno.txt.gz) | canonical_prediction_and_criteria_match |
| InterVar_2018 | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38_fixed.hg38_multianno.txt.gz](../data/tool_results/InterVar_2018/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38_fixed.hg38_multianno.txt.gz) | canonical_prediction_and_criteria_match |
| InterVar_2018 | clingen_28012026 | [clingen_28012026_hg38.hg38_multianno.txt.gz](../data/tool_results/InterVar_2018/clingen_28012026/clingen_28012026_hg38.hg38_multianno.txt.gz) | canonical_prediction_and_criteria_match |
| InterVar_2018 | foxl2 | [foxl2_hg38_intervar_20180118.txt.gz](../data/tool_results/InterVar_2018/foxl2/foxl2_hg38_intervar_20180118.txt.gz) | canonical_prediction_and_criteria_match |
| InterVar_2025 | HGMD_Clinvar_Cancer | [HGMD_Clinvar_Cancer_hg38_fixed.hg38_multianno.txt.gz](../data/tool_results/InterVar_2025/HGMD_Clinvar_Cancer/HGMD_Clinvar_Cancer_hg38_fixed.hg38_multianno.txt.gz) | canonical_prediction_and_criteria_match |
| InterVar_2025 | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg38_fixed.hg38_multianno.txt.gz](../data/tool_results/InterVar_2025/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg38_fixed.hg38_multianno.txt.gz) | canonical_prediction_and_criteria_match |
| InterVar_2025 | clingen_28012026 | [clingen_28012026_hg38.hg38_multianno.txt.gz](../data/tool_results/InterVar_2025/clingen_28012026/clingen_28012026_hg38.hg38_multianno.txt.gz) | canonical_prediction_and_criteria_match |
| InterVar_2025 | foxl2 | [foxl2_hg38_intervar_20250721.txt.gz](../data/tool_results/InterVar_2025/foxl2/foxl2_hg38_intervar_20250721.txt.gz) | canonical_prediction_and_criteria_match |
| TAPES | HGMD_Clinvar_Cancer | [tapes_HGMD_Clinvar_Cancer_hg38.csv.gz](../data/tool_results/TAPES/HGMD_Clinvar_Cancer/tapes_HGMD_Clinvar_Cancer_hg38.csv.gz) | canonical_prediction_and_criteria_match |
| TAPES | HGMD_Clinvar_HL | [tapes_HGMD_Clinvar_HL_hg38.csv.gz](../data/tool_results/TAPES/HGMD_Clinvar_HL/tapes_HGMD_Clinvar_HL_hg38.csv.gz) | canonical_prediction_and_criteria_match |
| TAPES | clingen_28012026 | [tapes_clingen_28012026_hg38_fixed.csv.gz](../data/tool_results/TAPES/clingen_28012026/tapes_clingen_28012026_hg38_fixed.csv.gz) | canonical_prediction_and_criteria_match |
| TAPES | foxl2 | [foxl2_hg38_tapes.csv.gz](../data/tool_results/TAPES/foxl2/foxl2_hg38_tapes.csv.gz) | canonical_prediction_and_criteria_match |
| VIP-HL | HGMD_Clinvar_HL | [HGMD_Clinvar_HL_hg19_viphl_results.tsv.gz](../data/tool_results/VIP-HL/HGMD_Clinvar_HL/HGMD_Clinvar_HL_hg19_viphl_results.tsv.gz) | canonical_prediction_and_criteria_match |
| VIP-HL | clingen_28012026 | [clingen_hearing_loss_28012026_hg19_viphl_results.tsv.gz](../data/tool_results/VIP-HL/clingen_28012026/clingen_hearing_loss_28012026_hg19_viphl_results.tsv.gz) | canonical_prediction_and_criteria_match |
| VIP-HL | foxl2 | [foxl2_hg19_viphl.tsv.gz](../data/tool_results/VIP-HL/foxl2/foxl2_hg19_viphl.tsv.gz) | criteria_match_eight_unknown_vs_missing_classifications |
