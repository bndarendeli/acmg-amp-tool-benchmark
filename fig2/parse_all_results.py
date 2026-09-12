#!/usr/bin/env python3
"""
Parse all tool results from Results folder and merge with ground truth.

This script:
1. Loads master ground truth (deduplicated)
2. Parses all tool results for all datasets
3. Merges tool results with ground truth
4. Saves merged results for: all variants, per dataset
"""

import pandas as pd
import sys
from pathlib import Path

# Add parsers to path
sys.path.insert(0, '/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_2/parsers')
from parser_factory import ParserFactory

# Directories
RESULTS_DIR = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/Results')
GT_DIR = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_2/ground_truth')
OUTPUT_DIR = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_2/results')

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

# Dataset files to process (all 5 main datasets)
DATASETS = [
    'clingen_28012026',
    'clingen_cancerpredisposition_28012026',
    'clingen_hearing_loss_28012026',
    'HGMD_Clinvar_Cancer',
    'HGMD_Clinvar_HL'
]


def make_variant_key(chrom, pos, ref, alt):
    """Create standardized variant key"""
    return f"{chrom}-{pos}-{ref}-{alt}"


def parse_tool_results(tool_folder, tool_name, dataset_name):
    """Parse results for a specific tool and dataset"""
    tool_dir = RESULTS_DIR / tool_folder
    
    # Tool-specific file patterns
    patterns = []
    
    if tool_name in ['InterVar_2018', 'InterVar_2025']:
        # InterVar: clingen_28012026_hg38.hg38_multianno.txt
        patterns = [
            f"{dataset_name}_hg38.hg38_multianno.txt",
            f"{dataset_name.replace('clingen_', 'clingen_hg38_')}.hg38_multianno.txt",
            f"{dataset_name}_hg38_fixed.hg38_multianno.txt"
        ]
    elif tool_name == 'Genebe':
        # Genebe: clingen_28012026_hg38_genebe.vcf
        patterns = [
            f"{dataset_name}_hg38_genebe.vcf",
            f"{dataset_name.replace('clingen_', 'clingen_')}_hg38_genebe.vcf"
        ]
    elif tool_name == 'Franklin':
        # Franklin: merged_all_clingen_28012026_hg19_last_single_snp_variants.csv
        patterns = [
            f"merged_all_{dataset_name}_hg19_last_single_snp_variants.csv",
            f"merged_{dataset_name}_hg19_last_single_snp_variants.csv",
            f"merged_variants-{dataset_name}_hg19_last.csv",
            f"merged_variants-{dataset_name.replace('clingen_hearing_loss_28012026', 'clingen_hl')}.csv"
        ]
    elif tool_name == 'TAPES':
        # TAPES: tapes_clingen_28012026_hg38_fixed.csv
        patterns = [
            f"tapes_{dataset_name}_hg38_fixed.csv",
            f"tapes_{dataset_name}_hg38.csv"
        ]
    elif tool_name == 'DiabloACMG':
        # DiabloACMG: clingen_hg38_28012026_annotated.tsv
        patterns = [
            f"clingen_hg38_{dataset_name.replace('clingen_', '')}_annotated.tsv",
            f"{dataset_name.replace('clingen_', 'clingen_hg38_')}_annotated.tsv",
            f"HGMD_hg38_{dataset_name.replace('HGMD_Clinvar_', 'Clinvar_')}_annotated.tsv"
        ]
    elif tool_name == 'VIP-HL':
        # VIP-HL: HGMD_Clinvar_HL_hg19_viphl_results.tsv, clingen_hearing_loss_28012026_hg19_viphl_results.tsv
        patterns = [
            f"{dataset_name}_hg19_viphl_results.tsv",
            f"{dataset_name.replace('clingen_', '')}_hg19_viphl_results.tsv",
            f"{dataset_name}_hg19_viphl.tsv"
        ]
    elif tool_name == 'CancerSIGVAR':
        # CancerSIGVAR: HGMD_Clinvar_Cancer_hg19_CancerSIGVAR_input_formatted_Cancer_SIGVAR_Results.tsv
        patterns = [
            f"{dataset_name}_hg19_CancerSIGVAR_input_formatted_Cancer_SIGVAR_Results.tsv",
            f"{dataset_name.replace('clingen_cancerpredisposition_28012026', 'clingen_cancerpredisposition_28012026')}_hg19_CancerSIGVAR_input_formatted_Cancer_SIGVAR_Results.tsv"
        ]
    elif tool_name == 'CPSR':
        # CPSR: HGMD_Clinvar_Cancer_hg38_grch38.cpsr.grch38.classification.tsv
        patterns = [
            f"{dataset_name}_hg38_grch38.cpsr.grch38.classification.tsv",
            f"{dataset_name.replace('clingen_cancerpredisposition_28012026', 'clingen_cancerpred')}_grch38.cpsr.grch38.classification.tsv",
            f"{dataset_name.replace('clingen_hearing_loss_28012026', 'clingen_hl')}_grch38.cpsr.grch38.classification.tsv"
        ]
    else:
        # Default patterns for other tools
        patterns = [
            f"{dataset_name}_hg38.tsv",
            f"{dataset_name}_hg38.vcf",
            f"{dataset_name}_hg38*.tsv"
        ]
    
    result_file = None
    for pattern in patterns:
        files = list(tool_dir.glob(pattern))
        if files:
            result_file = files[0]
            break
    
    if not result_file:
        return None
    
    # Get appropriate parser
    parser = ParserFactory.get_parser(tool_name, str(result_file))
    if not parser:
        print(f"  Warning: No parser for {tool_name}")
        return None
    
    # Parse the file
    try:
        df = parser.parse(str(result_file))
        if df is not None and len(df) > 0:
            # Deduplicate
            df = df.drop_duplicates(subset='Variant_Key', keep='first')
            return df
    except Exception as e:
        print(f"  Error parsing {tool_name} for {dataset_name}: {e}")
    
    return None


def main():
    print("="*80)
    print("PARSING ALL RESULTS")
    print("="*80)
    
    # Load master ground truth (5 datasets kept separate)
    print("\nLoading master ground truth...")
    master_gt = pd.read_excel(GT_DIR / 'master_ground_truth_5datasets.xlsx')
    
    # Create hg38 variant key (default)
    master_gt['Variant_Key'] = (master_gt['Chr'].astype(str) + '-' + 
                                 master_gt['Pos'].astype(str) + '-' + 
                                 master_gt['Ref'].astype(str) + '-' + 
                                 master_gt['Alt'].astype(str))
    
    # Create hg19 variant key for Franklin and VIP-HL
    master_gt['Variant_Key_hg19'] = (master_gt['Chr_hg19'].astype(str) + '-' + 
                                      master_gt['Pos_hg19'].astype(str) + '-' + 
                                      master_gt['Ref'].astype(str) + '-' + 
                                      master_gt['Alt'].astype(str))
    
    print(f"  Loaded {len(master_gt)} unique variants")
    print(f"  hg38 coordinates: {master_gt['Variant_Key'].notna().sum()}")
    print(f"  hg19 coordinates: {master_gt['Variant_Key_hg19'].notna().sum()}")
    
    # Initialize merged dataframe with ground truth
    merged_all = master_gt[['Variant_Key', 'Variant_Key_hg19', 'Chr', 'Pos', 
                            'Chr_hg19', 'Pos_hg19', 'Ref', 'Alt', 
                            'CLASSIFICATION', 'MET_CODES', 'Dataset']].copy()
    merged_all = merged_all.rename(columns={
        'CLASSIFICATION': 'Ground_Truth_Classification',
        'MET_CODES': 'Ground_Truth_ACMG'
    })
    
    # Parse all tools for all datasets
    print("\nParsing tool results...")
    for tool_folder, tool_name in TOOL_MAPPING.items():
        print(f"\n{tool_name}:")
        
        # Collect results from all datasets for this tool
        tool_results = []
        for dataset in DATASETS:
            df = parse_tool_results(tool_folder, tool_name, dataset)
            if df is not None:
                print(f"  {dataset}: {len(df)} variants")
                tool_results.append(df)
        
        if tool_results:
            # Combine all datasets for this tool
            combined = pd.concat(tool_results, ignore_index=True)
            # Deduplicate across datasets
            combined = combined.drop_duplicates(subset='Variant_Key', keep='first')
            print(f"  Total unique: {len(combined)} variants")
            
            # Determine which variant key to use for merging
            # Franklin, VIP-HL, and CancerSIGVAR use hg19 coordinates
            if tool_name in ['Franklin', 'VIP-HL', 'CancerSIGVAR']:
                merge_key_gt = 'Variant_Key_hg19'
                merge_key_tool = 'Variant_Key'
                print(f"  Using hg19 coordinates for {tool_name}")
            else:
                merge_key_gt = 'Variant_Key'
                merge_key_tool = 'Variant_Key'
            
            # Merge with master
            merged_all = merged_all.merge(
                combined[['Variant_Key', 'Classification', 'ACMG_Criteria']],
                left_on=merge_key_gt,
                right_on=merge_key_tool,
                how='left',
                suffixes=('', f'_{tool_name}')
            )
            
            # Drop the extra Variant_Key column from tool if it was created
            if f'Variant_Key_{tool_name}' in merged_all.columns:
                merged_all = merged_all.drop(columns=[f'Variant_Key_{tool_name}'])
            
            merged_all = merged_all.rename(columns={
                'Classification': f'{tool_name}_Classification',
                'ACMG_Criteria': f'{tool_name}_ACMG_Criteria'
            })
    
    # Save master merged results
    print("\n\nSaving results...")
    output_file = OUTPUT_DIR / 'all_variants' / 'merged_results.tsv'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    merged_all.to_csv(output_file, sep='\t', index=False)
    print(f"✓ Saved master merged results: {output_file}")
    print(f"  Total variants: {len(merged_all)}")
    print(f"  Total columns: {len(merged_all.columns)}")
    
    # Create per-dataset merged results
    print("\nCreating per-dataset merged results...")
    for dataset in DATASETS:
        dataset_merged = merged_all[merged_all['Dataset'] == dataset].copy()
        if len(dataset_merged) > 0:
            output_file = OUTPUT_DIR / 'per_dataset' / f'{dataset}_merged_results.tsv'
            output_file.parent.mkdir(parents=True, exist_ok=True)
            dataset_merged.to_csv(output_file, sep='\t', index=False)
            print(f"  {dataset}: {len(dataset_merged)} variants")
    
    print("\n" + "="*80)
    print("PARSING COMPLETE")
    print("="*80)


if __name__ == '__main__':
    main()
