#!/usr/bin/env python3
"""
Recalculate Call Rates Correctly

Fixes the call rate calculation by counting actual classified variants
per dataset, not using the Total_Variants from metrics summary which
includes all datasets combined.

Correct Call Rate = (Non-null classifications in dataset) / (Total variants in dataset)
"""

import pandas as pd
from pathlib import Path as FilePath


def recalculate_call_rates_for_dataset(data_file, dataset_id, dataset_name, total_variants, tools):
    """Calculate correct call rates for a specific dataset"""
    
    print(f"\n  Processing {dataset_name}...")
    
    # Load data
    df = pd.read_csv(data_file, sep='\t', low_memory=False)
    dataset_df = df[df['Dataset'] == dataset_id]
    
    if len(dataset_df) == 0:
        print(f"    ⚠ No data found for {dataset_id}")
        return None
    
    print(f"    Ground truth variants: {len(dataset_df):,}")
    
    results = []
    
    for tool in tools:
        # Find classification column for this tool
        classification_col = f'{tool}_Classification'
        
        if classification_col not in dataset_df.columns:
            continue
        
        # Count actual classifications (excluding 'Unknown', NULL, and empty)
        # Variant Resolution Rate = Only P, LP, VUS, LB, B count as resolved
        classified_mask = (
            dataset_df[classification_col].notna() & 
            (dataset_df[classification_col] != '') & 
            (dataset_df[classification_col] != 'Unknown')
        )
        classified = classified_mask.sum()
        
        call_rate = classified / total_variants * 100
        
        results.append({
            'Tool': tool,
            'Dataset': dataset_name,
            'Total_Variants': total_variants,
            'Classified': classified,
            'Call_Rate': call_rate
        })
        
        print(f"    {tool:<20} {classified:>6}/{total_variants:<6} = {call_rate:>6.2f}%")
    
    return pd.DataFrame(results)


def main():
    """Main execution"""
    
    print("=" * 80)
    print("RECALCULATING CALL RATES CORRECTLY")
    print("=" * 80)
    
    base_dir = FilePath(__file__).parent.parent
    data_file = base_dir / 'data' / 'merged_results.tsv'
    output_dir = base_dir / 'results' / 'corrected_metrics'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # All tools
    tools = [
        'InterVar_2018', 'InterVar_2025', 'BIAS', 'CharGer_Local', 'CharGer_Online',
        'DiabloACMG', 'Exomiser', 'Genebe', 'TAPES', 'Franklin', 'AutoGVP',
        'VIP-HL', 'CancerSIGVAR', 'CPSR'
    ]
    
    # Datasets
    datasets = [
        ('clingen_28012026', 'ClinGen', 11409),
        ('HGMD_Clinvar_HL', 'HGMD+ClinVar HL', 1948),
        ('HGMD_Clinvar_Cancer', 'HGMD+ClinVar Cancer', 1788)
    ]
    
    print("\n1. Calculating correct call rates from actual data...")
    
    all_results = []
    
    for dataset_id, dataset_name, total in datasets:
        result_df = recalculate_call_rates_for_dataset(
            data_file, dataset_id, dataset_name, total, tools
        )
        if result_df is not None:
            all_results.append(result_df)
    
    # Combine all results
    combined = pd.concat(all_results, ignore_index=True)
    
    # Add FOXL2 from separate analysis
    print(f"\n  Processing FOXL2 (from separate analysis)...")
    foxl2_dir = base_dir / 'results' / 'foxl2'
    foxl2_metrics = foxl2_dir / 'classification_metrics_summary.csv'
    
    if foxl2_metrics.exists():
        foxl2_df = pd.read_csv(foxl2_metrics)
        foxl2_df['Dataset'] = 'FOXL2'
        foxl2_df['Call_Rate'] = foxl2_df['Total_Variants'] / 288 * 100
        foxl2_df = foxl2_df.rename(columns={'Total_Variants': 'Classified'})
        foxl2_df['Total_Variants'] = 288
        
        foxl2_results = foxl2_df[['Tool', 'Dataset', 'Total_Variants', 'Classified', 'Call_Rate']]
        combined = pd.concat([combined, foxl2_results], ignore_index=True)
        
        print(f"    Added FOXL2 data for {len(foxl2_results)} tools")
    
    # Save corrected call rates
    output_file = output_dir / 'corrected_call_rates.csv'
    combined.to_csv(output_file, index=False)
    
    print(f"\n2. Saved corrected call rates to: {output_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY - Corrected Call Rates")
    print("=" * 80)
    
    for dataset in combined['Dataset'].unique():
        dataset_df = combined[combined['Dataset'] == dataset]
        print(f"\n{dataset}:")
        print(f"  Mean Call Rate: {dataset_df['Call_Rate'].mean():.2f}%")
        print(f"  Tools with 100% coverage: {len(dataset_df[dataset_df['Call_Rate'] == 100])}")
        print(f"  Tools with >100% coverage: {len(dataset_df[dataset_df['Call_Rate'] > 100])}")
        print(f"  Tools with <50% coverage: {len(dataset_df[dataset_df['Call_Rate'] < 50])}")
    
    print("\n" + "=" * 80)
    print("✓ Call rates recalculated correctly!")
    print("=" * 80)
    
    return combined


if __name__ == '__main__':
    main()
