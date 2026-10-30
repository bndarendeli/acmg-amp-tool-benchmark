import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

import seaborn as sns

from pathlib import Path

from sklearn.metrics import accuracy_score, precision_recall_fscore_support, confusion_matrix

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

def load_foxl2_data():
    """Load FOXL2 merged results"""
    base_dir = Path(__file__).resolve().parents[2]
    foxl2_file = base_dir / 'data/processed/foxl2_merged_results.tsv'
    
    if not foxl2_file.exists():
        print(f"✗ FOXL2 data file not found: {foxl2_file}")
        return None
    
    df = pd.read_csv(foxl2_file, sep='\t')
    print(f"✓ Loaded {len(df)} FOXL2 variants")
    
    return df

def calculate_jaccard_per_criterion(df):
    """Calculate Jaccard similarity per tool per criterion"""
    results = []
    
    for tool in TOOLS:
        criteria_col = f'{tool}_ACMG_Criteria'
        if criteria_col not in df.columns:
            continue
        
        for criterion in ALL_CRITERIA:
            gt_has_criterion = []
            tool_has_criterion = []
            
            for idx, row in df.iterrows():
                gt_criteria = parse_criteria_string(row.get('Ground_Truth_ACMG', ''))
                tool_criteria = parse_criteria_string(row.get(criteria_col, ''))
                
                gt_has = criterion in gt_criteria
                tool_has = criterion in tool_criteria
                
                gt_has_criterion.append(gt_has)
                tool_has_criterion.append(tool_has)
            
            gt_set = set([i for i, x in enumerate(gt_has_criterion) if x])
            tool_set = set([i for i, x in enumerate(tool_has_criterion) if x])
            
            intersection = len(gt_set.intersection(tool_set))
            union = len(gt_set.union(tool_set))
            jaccard = intersection / union if union > 0 else 0.0
            
            results.append({
                'Tool': tool,
                'Criterion': criterion,
                'Jaccard_Similarity': jaccard,
                'Intersection': intersection,
                'Union': union
            })
    
    return pd.DataFrame(results)

def main():
    root = Path(__file__).resolve().parents[2]
    df = pd.read_csv(root / 'data/processed/foxl2_merged_results.tsv', sep='\t', low_memory=False)
    result = calculate_jaccard_per_criterion(df)
    result.to_csv(Path(__file__).resolve().parent / 'jaccard_per_criterion_foxl2.csv', index=False)

if __name__ == '__main__':
    main()
