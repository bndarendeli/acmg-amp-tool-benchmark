#!/usr/bin/env python3
"""
Build Figure 5 Master Table: Evidence-Level Concordance Analysis

Generates master table with one row per dataset × variant × tool.
Calculates classification concordance, evidence concordance, and Jaccard similarity.

Only includes ClinGen and FOXL2 datasets (datasets with reference ACMG criteria).
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Import existing criterion parser
sys.path.insert(0, str(Path(__file__).parent.parent / 'fig1'))
from utils_criteria_parser import parse_criteria_string

# Configuration
BASE_DIR = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep')
DATA_FILE = BASE_DIR / 'analysis_last' / 'data' / 'merged_results.tsv'
FOXL2_FILE = BASE_DIR / 'analysis_last' / 'foxl2_merged_results.tsv'
OUTPUT_DIR = Path(__file__).parent / 'outputs'
OUTPUT_FILE = OUTPUT_DIR / 'fig5_master_variant_tool.csv'

# Only datasets with reference ACMG criteria
INCLUDED_DATASETS = ['clingen_28012026', 'foxl2']

# All tools
TOOLS = [
    'InterVar_2018', 'InterVar_2025', 'BIAS', 'CharGer_Local', 'CharGer_Online',
    'DiabloACMG', 'Exomiser', 'Genebe', 'TAPES', 'Franklin', 'AutoGVP', 'VIP-HL',
    'CancerSIGVAR', 'CPSR'
]

# Classification normalization mapping
# Handles both full names (ClinGen) and abbreviations (FOXL2)
CLASS_NORMALIZATION = {
    # Full names (ClinGen format)
    'Pathogenic': 'P',
    'Likely_Pathogenic': 'LP',
    'Uncertain_Significance': 'VUS',
    'Likely_Benign': 'LB',
    'Benign': 'B',
    # Abbreviations (FOXL2 format)
    'P': 'P',
    'LP': 'LP',
    'VUS': 'VUS',
    'LB': 'LB',
    'B': 'B',
    # Excluded values
    'Unknown': None,
    'Not_Provided': None,
    '': None,
    None: None
}


def normalize_classification(classification):
    """Normalize classification to five-class system"""
    if pd.isna(classification):
        return None
    
    classification = str(classification).strip()
    return CLASS_NORMALIZATION.get(classification, None)


def calculate_jaccard(set_a, set_b, a_available=True, b_available=True):
    """
    Calculate Jaccard similarity between two criterion sets.
    
    Rules:
    - Both sets non-empty: standard Jaccard
    - Reference empty, tool non-empty: 0
    - Reference non-empty, tool empty: 0
    - Both truly observed and empty: 1
    - Either side missing/unavailable: NA
    
    Args:
        set_a: First criterion set
        set_b: Second criterion set
        a_available: Whether set_a data is available (not missing)
        b_available: Whether set_b data is available (not missing)
    
    Returns:
        Jaccard similarity [0, 1] or None (NA)
    """
    # If either side is missing/unavailable, return NA
    if not a_available or not b_available:
        return None
    
    # Both available - calculate Jaccard
    if len(set_a) == 0 and len(set_b) == 0:
        # Both truly observed and empty
        return 1.0
    
    union = set_a.union(set_b)
    if len(union) == 0:
        return 1.0
    
    intersection = set_a.intersection(set_b)
    return len(intersection) / len(union)


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


def build_master_table(df):
    """
    Build master table with one row per dataset × variant × tool.
    
    Returns:
        DataFrame with evidence-level concordance metrics
    """
    print("\n" + "="*80)
    print("BUILDING FIGURE 4 MASTER TABLE")
    print("="*80)
    
    master_rows = []
    
    # Filter to included datasets only
    df_filtered = df[df['Dataset'].isin(INCLUDED_DATASETS)].copy()
    print(f"\nFiltered to {len(df_filtered)} variants from {INCLUDED_DATASETS}")
    
    for dataset in INCLUDED_DATASETS:
        dataset_df = df_filtered[df_filtered['Dataset'] == dataset]
        print(f"\nProcessing {dataset}: {len(dataset_df)} variants")
        
        for tool in TOOLS:
            tool_class_col = f'{tool}_Classification'
            tool_criteria_col = f'{tool}_ACMG_Criteria'
            
            # Check if tool columns exist
            if tool_class_col not in df.columns:
                print(f"  ⚠ {tool}: Classification column not found")
                continue
            
            tool_has_criteria = tool_criteria_col in df.columns
            if not tool_has_criteria:
                print(f"  ⚠ {tool}: Criteria column not found")
            
            n_processed = 0
            n_class_available = 0
            n_criteria_available = 0
            
            for idx, row in dataset_df.iterrows():
                # Get variant identifiers
                variant_id = row.get('Variant_Key', '')
                chr_val = row.get('Chr', '')
                pos_val = row.get('Pos', '')
                ref_val = row.get('Ref', '')
                alt_val = row.get('Alt', '')
                
                # Get and normalize classifications
                ref_class_raw = row.get('Ground_Truth_Classification')
                tool_class_raw = row.get(tool_class_col)
                
                ref_class = normalize_classification(ref_class_raw)
                tool_class = normalize_classification(tool_class_raw)
                
                # Skip if either classification is missing or Unknown
                if ref_class is None or tool_class is None:
                    continue
                
                n_class_available += 1
                
                # Get reference criteria
                ref_criteria_raw = row.get('Ground_Truth_ACMG', '')
                ref_criteria_available = is_criteria_available(ref_criteria_raw)
                ref_criteria_set = parse_criteria_string(ref_criteria_raw) if ref_criteria_available else set()
                
                # Get tool criteria
                tool_criteria_raw = row.get(tool_criteria_col, '') if tool_has_criteria else ''
                tool_criteria_available = is_criteria_available(tool_criteria_raw) if tool_has_criteria else False
                tool_criteria_set = parse_criteria_string(tool_criteria_raw) if tool_criteria_available else set()
                
                if ref_criteria_available and tool_criteria_available:
                    n_criteria_available += 1
                
                # Calculate metrics
                classification_match = (ref_class == tool_class)
                
                # Evidence concordance only if both criteria are available
                if ref_criteria_available and tool_criteria_available:
                    exact_evidence_match = (ref_criteria_set == tool_criteria_set)
                    shared = ref_criteria_set.intersection(tool_criteria_set)
                    missing = ref_criteria_set.difference(tool_criteria_set)
                    additional = tool_criteria_set.difference(ref_criteria_set)
                    
                    jaccard = calculate_jaccard(
                        ref_criteria_set, tool_criteria_set,
                        ref_criteria_available, tool_criteria_available
                    )
                else:
                    exact_evidence_match = None
                    shared = set()
                    missing = set()
                    additional = set()
                    jaccard = None
                
                # Build row
                master_rows.append({
                    'dataset': dataset,
                    'variant_id': variant_id,
                    'chr': chr_val,
                    'pos': pos_val,
                    'ref': ref_val,
                    'alt': alt_val,
                    'tool': tool,
                    'reference_class': ref_class,
                    'tool_class': tool_class,
                    'reference_criteria': ','.join(sorted(ref_criteria_set)) if ref_criteria_available else None,
                    'tool_criteria': ','.join(sorted(tool_criteria_set)) if tool_criteria_available else None,
                    'classification_match': classification_match,
                    'exact_evidence_match': exact_evidence_match,
                    'criterion_jaccard': jaccard,
                    'shared_criteria': ','.join(sorted(shared)) if (ref_criteria_available and tool_criteria_available) else None,
                    'missing_criteria': ','.join(sorted(missing)) if (ref_criteria_available and tool_criteria_available) else None,
                    'additional_criteria': ','.join(sorted(additional)) if (ref_criteria_available and tool_criteria_available) else None,
                    'n_reference_criteria': len(ref_criteria_set) if ref_criteria_available else None,
                    'n_tool_criteria': len(tool_criteria_set) if tool_criteria_available else None,
                    'n_shared': len(shared) if (ref_criteria_available and tool_criteria_available) else None,
                    'n_missing': len(missing) if (ref_criteria_available and tool_criteria_available) else None,
                    'n_additional': len(additional) if (ref_criteria_available and tool_criteria_available) else None,
                })
                
                n_processed += 1
            
            if n_processed > 0:
                print(f"  ✓ {tool}: {n_processed} variants, {n_class_available} with classifications, {n_criteria_available} with criteria")
    
    master_df = pd.DataFrame(master_rows)
    
    print(f"\n{'='*80}")
    print(f"MASTER TABLE CREATED")
    print(f"{'='*80}")
    print(f"Total rows: {len(master_df):,}")
    print(f"Datasets: {master_df['dataset'].nunique()}")
    print(f"Unique variants: {master_df['variant_id'].nunique()}")
    print(f"Tools: {master_df['tool'].nunique()}")
    print(f"Rows with classification data: {master_df['classification_match'].notna().sum():,}")
    print(f"Rows with evidence data: {master_df['exact_evidence_match'].notna().sum():,}")
    
    return master_df


def load_all_data():
    """
    Load data from multiple sources and combine.
    
    Returns:
        Combined DataFrame with all datasets
    """
    print("\n" + "="*80)
    print("LOADING DATA FROM MULTIPLE SOURCES")
    print("="*80)
    
    all_dfs = []
    
    # Load main merged_results.tsv (ClinGen and other datasets)
    print(f"\n1. Loading main data file: {DATA_FILE}")
    if DATA_FILE.exists():
        df_main = pd.read_csv(DATA_FILE, sep='\t', low_memory=False)
        print(f"   ✓ Loaded {len(df_main):,} variants")
        if 'Dataset' in df_main.columns:
            print(f"   Datasets: {df_main['Dataset'].unique().tolist()}")
        all_dfs.append(df_main)
    else:
        print(f"   ✗ File not found: {DATA_FILE}")
    
    # Load FOXL2 data
    print(f"\n2. Loading FOXL2 data: {FOXL2_FILE}")
    if FOXL2_FILE.exists():
        df_foxl2 = pd.read_csv(FOXL2_FILE, sep='\t', low_memory=False)
        print(f"   ✓ Loaded {len(df_foxl2):,} variants")
        
        # Add Dataset column if not present
        if 'Dataset' not in df_foxl2.columns:
            df_foxl2['Dataset'] = 'foxl2'
            print(f"   Added 'Dataset' column with value 'foxl2'")
        
        all_dfs.append(df_foxl2)
    else:
        print(f"   ✗ File not found: {FOXL2_FILE}")
        print(f"   Searching for alternative FOXL2 locations...")
        
        # Try alternative locations
        alternative_paths = [
            BASE_DIR / 'datasets' / 'foxl2' / 'analysis' / 'results' / 'foxl2_merged_results.tsv',
            BASE_DIR / 'analysis_2' / 'foxl2_merged_results.tsv',
            BASE_DIR / 'foxl2_merged_results.tsv',
        ]
        
        for alt_path in alternative_paths:
            if alt_path.exists():
                print(f"   ✓ Found FOXL2 data at: {alt_path}")
                df_foxl2 = pd.read_csv(alt_path, sep='\t', low_memory=False)
                print(f"   ✓ Loaded {len(df_foxl2):,} variants")
                
                if 'Dataset' not in df_foxl2.columns:
                    df_foxl2['Dataset'] = 'foxl2'
                    print(f"   Added 'Dataset' column with value 'foxl2'")
                
                all_dfs.append(df_foxl2)
                break
        else:
            print(f"   ⚠ FOXL2 data not found in any location")
    
    if not all_dfs:
        print("\n✗ ERROR: No data files found!")
        return None
    
    # Combine all dataframes
    print(f"\n{'='*80}")
    print("COMBINING DATA")
    print(f"{'='*80}")
    
    df_combined = pd.concat(all_dfs, ignore_index=True)
    print(f"✓ Combined {len(df_combined):,} total variants")
    
    if 'Dataset' in df_combined.columns:
        dataset_counts = df_combined['Dataset'].value_counts()
        print(f"\nDataset breakdown:")
        for dataset, count in dataset_counts.items():
            print(f"  {dataset}: {count:,} variants")
    
    return df_combined


def main():
    """Main execution"""
    print("="*80)
    print("FIGURE 4: EVIDENCE-LEVEL CONCORDANCE ANALYSIS")
    print("Master Table Builder")
    print("="*80)
    
    # Load data from multiple sources
    df = load_all_data()
    if df is None:
        return 1
    
    # Build master table
    master_df = build_master_table(df)
    
    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    master_df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✓ Master table saved to: {OUTPUT_FILE}")
    print(f"  File size: {OUTPUT_FILE.stat().st_size / 1024:.1f} KB")
    
    # Quick statistics
    print(f"\n{'='*80}")
    print("QUICK STATISTICS")
    print(f"{'='*80}")
    
    print("\nBy Dataset:")
    for dataset in master_df['dataset'].unique():
        dataset_df = master_df[master_df['dataset'] == dataset]
        print(f"  {dataset}:")
        print(f"    Rows: {len(dataset_df):,}")
        print(f"    Unique variants: {dataset_df['variant_id'].nunique()}")
        print(f"    With evidence data: {dataset_df['exact_evidence_match'].notna().sum():,}")
    
    print("\nClassification Concordance:")
    total_with_class = master_df['classification_match'].notna().sum()
    concordant = master_df['classification_match'].sum()
    print(f"  Total evaluated: {total_with_class:,}")
    print(f"  Concordant: {concordant:,} ({100*concordant/total_with_class:.1f}%)")
    
    print("\nEvidence Concordance (among class-concordant):")
    class_concordant = master_df[master_df['classification_match'] == True]
    with_evidence = class_concordant['exact_evidence_match'].notna().sum()
    evidence_concordant = class_concordant['exact_evidence_match'].sum()
    if with_evidence > 0:
        print(f"  Class-concordant with evidence data: {with_evidence:,}")
        print(f"  Evidence concordant: {evidence_concordant:,} ({100*evidence_concordant/with_evidence:.1f}%)")
    else:
        print(f"  No variants with both class concordance and evidence data")
    
    print(f"\n{'='*80}")
    print("✓ MASTER TABLE BUILD COMPLETE")
    print(f"{'='*80}\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
