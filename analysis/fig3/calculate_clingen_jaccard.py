import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

import seaborn as sns

from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'fig2'))
from utils_criteria_parser import parse_criteria_string

TOOLS = [
    'InterVar_2018', 'InterVar_2025', 'BIAS', 'CharGer_Local', 'CharGer_Online',
    'DiabloACMG', 'Exomiser', 'Genebe', 'TAPES', 'Franklin', 'AutoGVP', 'VIP-HL',
    'CancerSIGVAR', 'CPSR'
]

ALL_CRITERIA = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]

def is_criteria_available(criteria_str):
    """
    Determine if criteria data is available (not missing).
    
    Missing data indicators: 
    - None, NaN, empty string
    - 'Unknown', 'Not_Provided', 'NA', 'nan'
    - '.' (VCF-style missing value placeholder)
    
    CRITICAL: '.' in ClinGen data represents MISSING/UNAVAILABLE criterion
    information, NOT a genuinely observed empty criterion set.
    """
    if pd.isna(criteria_str):
        return False
    
    criteria_str = str(criteria_str).strip()
    if criteria_str in ['', '.', 'Unknown', 'Not_Provided', 'None', 'nan', 'NA']:
        return False
    
    return True

def calculate_jaccard_per_tool_criterion(df, dataset_name):
    """
    Calculate Jaccard similarity per tool per criterion.
    
    For each tool and criterion:
    - Find all variants where GT has the criterion
    - Find all variants where Tool has the criterion
    - Calculate Jaccard = |GT ∩ Tool| / |GT ∪ Tool|
    
    IMPORTANT: Excludes variants with missing/unavailable reference criteria
    (Ground_Truth_ACMG = '.') to avoid treating missing data as empty sets.
    
    Returns DataFrame with columns: Tool, Criterion, Dataset, Jaccard_Similarity, Intersection, Union
    """
    results = []
    
    for tool in TOOLS:
        criteria_col = f'{tool}_ACMG_Criteria'
        if criteria_col not in df.columns:
            continue
        
        for criterion in ALL_CRITERIA:
            gt_has_criterion = []
            tool_has_criterion = []
            
            for idx, row in df.iterrows():
                # CRITICAL FIX: Skip variants with missing reference criteria
                gt_acmg = row.get('Ground_Truth_ACMG', '')
                if not is_criteria_available(gt_acmg):
                    continue
                
                gt_criteria = parse_criteria_string(gt_acmg)
                tool_criteria = parse_criteria_string(row.get(criteria_col, ''))
                
                gt_has = criterion in gt_criteria
                tool_has = criterion in tool_criteria
                
                gt_has_criterion.append(gt_has)
                tool_has_criterion.append(tool_has)
            
            # Create sets of variant indices
            gt_set = set([i for i, x in enumerate(gt_has_criterion) if x])
            tool_set = set([i for i, x in enumerate(tool_has_criterion) if x])
            
            intersection = len(gt_set.intersection(tool_set))
            union = len(gt_set.union(tool_set))
            jaccard = intersection / union if union > 0 else 0.0
            
            results.append({
                'Tool': tool,
                'Criterion': criterion,
                'Dataset': dataset_name,
                'Jaccard_Similarity': jaccard,
                'Intersection': intersection,
                'Union': union
            })
    
    return pd.DataFrame(results)

def main():
    root = Path(__file__).resolve().parents[2]
    merged = pd.read_csv(root / 'data/processed/merged_results_noHGMD.tsv', sep='\t', low_memory=False)
    clingen = merged[merged['Dataset'].eq('clingen_28012026')]
    result = calculate_jaccard_per_tool_criterion(clingen, 'clingen_28012026')
    result.to_csv(Path(__file__).resolve().parent / 'jaccard_per_criterion_clingen_28012026.csv', index=False)

if __name__ == '__main__':
    main()
