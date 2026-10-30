#!/usr/bin/env python3
"""
Parse tool results for ClinGen and FOXL2 and merge with ground truth (no HGMD).

This script:
1. Loads ClinGen ground truth (master ground truth filtered to clingen_28012026)
2. Parses ClinGen tool results (incl. hearing loss / cancer predisposition subsets,
   which provide VIP-HL, CancerSIGVAR and CPSR results for ClinGen variants)
3. Loads FOXL2 ground truth and parses FOXL2 tool results
4. Saves a single merged table: merged_results_noHGMD.tsv

Usage:
    python3 parse_all_results.py [output_file]
"""

import pandas as pd
import sys
from pathlib import Path

# Add parsers to path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'parsers'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from unpack_tool_results import unpack_tool_results
from parser_factory import ParserFactory

# Directories
REPO_DIR = Path(__file__).resolve().parents[2]
BASE_DIR = unpack_tool_results(REPO_DIR)
RESULTS_DIR = BASE_DIR / 'Results'
GT_DIR = REPO_DIR / 'data/reference'
FOXL2_GT_FILE = REPO_DIR / 'data/reference/FOXL2_dataset.xlsx'
FOXL2_RESULTS_DIR = BASE_DIR / 'datasets' / 'foxl2' / 'results_fix'
OUTPUT_FILE = REPO_DIR / 'data/processed/merged_results_noHGMD.tsv'

# Tool mapping (tool folder name -> parser name)
TOOL_MAPPING = {
    'Intervar_20180118': 'InterVar_2018',
    'Intervar_20250721': 'InterVar_2025',
    'BIAS': 'BIAS',
    'CharGer_local': 'CharGer_Local',
    'CharGer_online': 'CharGer_Online',
    'DiabloACMG': 'DiabloACMG',
    'Exomiser': 'Exomiser',
    'Genebe': 'Genebe',
    'TAPES': 'TAPES',
    'Franklin': 'Franklin',
    'AutoGVP': 'AutoGVP',
    'VIP-HL': 'VIP-HL',
    'Cancer SIGVAR': 'CancerSIGVAR',
    'CPSR': 'CPSR'
}

# ClinGen ground truth label in master ground truth
CLINGEN_DATASET = 'clingen_28012026'

# ClinGen tool result file sets (subsets share variants with clingen_28012026)
CLINGEN_RESULT_SETS = [
    'clingen_28012026',
    'clingen_cancerpredisposition_28012026',
    'clingen_hearing_loss_28012026'
]

# FOXL2 tool result files (datasets/foxl2/results_fix)
FOXL2_TOOL_FILES = {
    'InterVar_2018': 'foxl2_hg38_intervar_20180118.txt',
    'InterVar_2025': 'foxl2_hg38_intervar_20250721.txt',
    'BIAS': 'foxl2_hg38_BIAS.tsv',
    'CharGer_Local': 'foxl2_hg38_CharGer_local_annotated_VEPv97.tsv',
    'CharGer_Online': 'foxl2_hg38_CharGer_online.tsv',
    'DiabloACMG': 'foxl2_hg38_DiabloACMG.tsv',
    'Exomiser': 'foxl2_hg38_exomiser.tsv',
    'Genebe': 'foxl2_hg38_genebe.vcf',
    'TAPES': 'foxl2_hg38_tapes.csv',
    'Franklin': 'foxl2_hg19_franklin.csv',
    'AutoGVP': 'foxl2_hg38_autogvp.tsv',
    'VIP-HL': 'foxl2_hg19_viphl.tsv'
}

# FOXL2 ground truth labels -> standard labels used in merged results
FOXL2_GT_MAP = {
    'P': 'Pathogenic',
    'LP': 'Likely_Pathogenic',
    'VUS': 'Uncertain_Significance',
    'LB': 'Likely_Benign',
    'B': 'Benign'
}

# Tools whose results use hg19 coordinates
HG19_TOOLS = ['Franklin', 'VIP-HL', 'CancerSIGVAR']

BASE_COLUMNS = ['Variant_Key', 'Variant_Key_hg19', 'Chr', 'Pos', 'Chr_hg19', 'Pos_hg19',
                'Ref', 'Alt', 'Ground_Truth_Classification', 'Ground_Truth_ACMG', 'Dataset']


def parse_file(tool_name, result_file):
    """Parse a single tool result file and deduplicate by Variant_Key"""
    parser = ParserFactory.get_parser(tool_name, str(result_file))
    if not parser:
        print(f"  Warning: No parser for {tool_name}")
        return None
    try:
        df = parser.parse(str(result_file))
        if df is not None and len(df) > 0:
            return df.drop_duplicates(subset='Variant_Key', keep='first')
    except Exception as e:
        print(f"  Error parsing {tool_name} ({result_file.name}): {e}")
    return None


def find_clingen_result_file(tool_folder, tool_name, dataset_name):
    """Locate the result file for a specific tool and ClinGen result set"""
    tool_dir = RESULTS_DIR / tool_folder

    if tool_name in ['InterVar_2018', 'InterVar_2025']:
        patterns = [
            f"{dataset_name}_hg38.hg38_multianno.txt",
            f"{dataset_name.replace('clingen_', 'clingen_hg38_')}.hg38_multianno.txt",
            f"{dataset_name}_hg38_fixed.hg38_multianno.txt"
        ]
    elif tool_name == 'Genebe':
        patterns = [f"{dataset_name}_hg38_genebe.vcf"]
    elif tool_name == 'Franklin':
        patterns = [
            f"merged_all_{dataset_name}_hg19_last_single_snp_variants.csv",
            f"merged_{dataset_name}_hg19_last_single_snp_variants.csv",
            f"merged_variants-{dataset_name}_hg19_last.csv",
            f"merged_variants-{dataset_name.replace('clingen_hearing_loss_28012026', 'clingen_hl')}.csv"
        ]
    elif tool_name == 'TAPES':
        patterns = [
            f"tapes_{dataset_name}_hg38_fixed.csv",
            f"tapes_{dataset_name}_hg38.csv"
        ]
    elif tool_name == 'DiabloACMG':
        patterns = [
            f"clingen_hg38_{dataset_name.replace('clingen_', '')}_annotated.tsv",
            f"{dataset_name.replace('clingen_', 'clingen_hg38_')}_annotated.tsv"
        ]
    elif tool_name == 'VIP-HL':
        patterns = [
            f"{dataset_name}_hg19_viphl_results.tsv",
            f"{dataset_name.replace('clingen_', '')}_hg19_viphl_results.tsv",
            f"{dataset_name}_hg19_viphl.tsv"
        ]
    elif tool_name == 'CancerSIGVAR':
        patterns = [f"{dataset_name}_hg19_CancerSIGVAR_input_formatted_Cancer_SIGVAR_Results.tsv"]
    elif tool_name == 'CPSR':
        patterns = [
            f"{dataset_name}_hg38_grch38.cpsr.grch38.classification.tsv",
            f"{dataset_name.replace('clingen_cancerpredisposition_28012026', 'clingen_cancerpred')}_grch38.cpsr.grch38.classification.tsv",
            f"{dataset_name.replace('clingen_hearing_loss_28012026', 'clingen_hl')}_grch38.cpsr.grch38.classification.tsv"
        ]
    else:
        patterns = [
            f"{dataset_name}_hg38.tsv",
            f"{dataset_name}_hg38.vcf",
            f"{dataset_name}_hg38*.tsv"
        ]

    for pattern in patterns:
        files = list(tool_dir.glob(pattern))
        if files:
            return files[0]
    return None


def merge_tool(merged, tool_df, tool_name):
    """Left-join one tool's classifications onto the ground truth table"""
    merge_key_gt = 'Variant_Key_hg19' if tool_name in HG19_TOOLS else 'Variant_Key'
    tool_data = tool_df[['Variant_Key', 'Classification', 'ACMG_Criteria']].rename(columns={
        'Variant_Key': '_tool_key',
        'Classification': f'{tool_name}_Classification',
        'ACMG_Criteria': f'{tool_name}_ACMG_Criteria'
    })
    merged = merged.merge(tool_data, left_on=merge_key_gt, right_on='_tool_key', how='left')
    return merged.drop(columns=['_tool_key'])


def build_clingen():
    """ClinGen ground truth merged with all tool results"""
    print("\n" + "=" * 80)
    print("CLINGEN")
    print("=" * 80)

    master_gt = pd.read_excel(GT_DIR / 'clingen_ground_truth.xlsx')
    gt = master_gt[master_gt['Dataset'] == CLINGEN_DATASET].copy()

    gt['Variant_Key'] = (gt['Chr'].astype(str) + '-' + gt['Pos'].astype(str) + '-' +
                         gt['Ref'].astype(str) + '-' + gt['Alt'].astype(str))
    gt['Variant_Key_hg19'] = (gt['Chr_hg19'].astype(str) + '-' + gt['Pos_hg19'].astype(str) + '-' +
                              gt['Ref'].astype(str) + '-' + gt['Alt'].astype(str))

    merged = gt[['Variant_Key', 'Variant_Key_hg19', 'Chr', 'Pos', 'Chr_hg19', 'Pos_hg19',
                 'Ref', 'Alt', 'CLASSIFICATION', 'MET_CODES', 'Dataset']].rename(columns={
        'CLASSIFICATION': 'Ground_Truth_Classification',
        'MET_CODES': 'Ground_Truth_ACMG'
    })
    print(f"  Ground truth variants: {len(merged)}")

    for tool_folder, tool_name in TOOL_MAPPING.items():
        print(f"\n{tool_name}:")
        tool_results = []
        for result_set in CLINGEN_RESULT_SETS:
            result_file = find_clingen_result_file(tool_folder, tool_name, result_set)
            if result_file is None:
                continue
            df = parse_file(tool_name, result_file)
            if df is not None:
                print(f"  {result_set}: {len(df)} variants")
                tool_results.append(df)

        if tool_results:
            combined = pd.concat(tool_results, ignore_index=True)
            combined = combined.drop_duplicates(subset='Variant_Key', keep='first')
            print(f"  Total unique: {len(combined)} variants")
            merged = merge_tool(merged, combined, tool_name)

    return merged


def build_foxl2():
    """FOXL2 ground truth merged with all tool results"""
    print("\n" + "=" * 80)
    print("FOXL2")
    print("=" * 80)

    raw = pd.read_excel(FOXL2_GT_FILE)

    def to_key(value):
        if pd.isna(value) or not str(value).strip():
            return ''
        parts = str(value).split('-')
        return '-'.join([parts[0].replace('chr', '')] + parts[1:]) if len(parts) == 4 else ''

    gt = pd.DataFrame({
        'Variant_Key': raw['hg38'].apply(to_key),
        'Variant_Key_hg19': raw['hg19'].apply(to_key),
        'Ground_Truth_Classification': raw['Classification'].map(FOXL2_GT_MAP).fillna(raw['Classification']),
        'Ground_Truth_ACMG': raw['ACMG-AMP classification'],
        'HGVS': raw['HGVS']
    })
    gt = gt[gt['Variant_Key'] != '']
    # Some variants have multiple HGVS notations but same genomic coordinates
    gt = gt.drop_duplicates(subset='Variant_Key', keep='first').reset_index(drop=True)

    hg38 = gt['Variant_Key'].str.split('-', n=3, expand=True)
    hg19 = gt['Variant_Key_hg19'].str.split('-', n=3, expand=True)
    gt['Chr'], gt['Pos'], gt['Ref'], gt['Alt'] = hg38[0], hg38[1].astype(int), hg38[2], hg38[3]
    gt['Chr_hg19'], gt['Pos_hg19'] = hg19[0], hg19[1].astype(int)
    gt['Dataset'] = 'FOXL2'

    merged = gt[BASE_COLUMNS + ['HGVS']]
    print(f"  Ground truth variants: {len(merged)}")

    for tool_name, filename in FOXL2_TOOL_FILES.items():
        result_file = FOXL2_RESULTS_DIR / filename
        if not result_file.exists():
            print(f"  Warning: {filename} not found, skipping {tool_name}")
            continue
        df = parse_file(tool_name, result_file)
        if df is not None:
            merged = merge_tool(merged, df, tool_name)
            print(f"  {tool_name}: {merged[f'{tool_name}_Classification'].notna().sum()} variants matched")

    return merged


def main():
    output_file = Path(sys.argv[1]) if len(sys.argv) > 1 else OUTPUT_FILE

    print("=" * 80)
    print("PARSING RESULTS (ClinGen + FOXL2, no HGMD)")
    print("=" * 80)

    clingen = build_clingen()
    foxl2 = build_foxl2()

    tool_columns = [c for t in TOOL_MAPPING.values()
                    for c in (f'{t}_Classification', f'{t}_ACMG_Criteria')]
    merged_all = pd.concat([clingen, foxl2], ignore_index=True)
    merged_all = merged_all.reindex(columns=BASE_COLUMNS + tool_columns + ['HGVS'])
    merged_all['Pos'] = merged_all['Pos'].astype('Int64')
    merged_all['Pos_hg19'] = merged_all['Pos_hg19'].astype('Int64')

    print("\n\nSaving results...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    merged_all.to_csv(output_file, sep='\t', index=False)
    print(f"✓ Saved merged results: {output_file}")
    print(f"  Total variants: {len(merged_all)}")
    print(f"  Per dataset: {merged_all['Dataset'].value_counts().to_dict()}")
    print(f"  Total columns: {len(merged_all.columns)}")

    print("\n" + "=" * 80)
    print("PARSING COMPLETE")
    print("=" * 80)


if __name__ == '__main__':
    main()
