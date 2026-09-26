#!/usr/bin/env python3
"""
Build Correlation Master Table for Criterion-Performance Analysis

Creates the authoritative ClinGen tool × metric master matrix for pairwise 
correlation analysis between criterion-level concordance and aggregate 
classification performance.

Output: 14 tools × 32 columns (1 Tool + 28 criteria + 3 aggregate metrics)
"""

import pandas as pd
import numpy as np
from pathlib import Path

print("=" * 80)
print("BUILDING CORRELATION MASTER TABLE")
print("=" * 80)

# Define paths
base_dir = Path(__file__).parent.parent
output_dir = Path(__file__).parent
output_dir.mkdir(parents=True, exist_ok=True)

# Source file paths
source_files = {
    'criterion_jaccard': base_dir / 'fig3' / 'jaccard_per_criterion_clingen_28012026.csv',
    'mean_jaccard': base_dir / 'fig5' / 'outputs' / 'fig4_capability_adjusted_summary.csv',
    'total_jaccard': base_dir / 'fig5' / 'outputs' / 'clingen_total_criterion_jaccard.csv',
    'accuracy': base_dir / 'fig4' / 'data' / 'reconstructed' / 'ClinGen_reconstructed_metrics.csv'
}

# ACMG/AMP criteria in specified order
CRITERIA_ORDER = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]

print("\n1. Loading source files...")
print("-" * 80)

# Verify all source files exist
for name, path in source_files.items():
    if not path.exists():
        raise FileNotFoundError(f"Source file not found: {path}")
    print(f"  ✓ {name}: {path}")

# Load criterion-specific Jaccard
print("\n2. Loading criterion-specific Jaccard...")
df_criterion = pd.read_csv(source_files['criterion_jaccard'])
print(f"  Loaded: {len(df_criterion)} rows")

# Validation: Check for exactly 392 rows (14 tools × 28 criteria)
expected_rows = 14 * 28
assert len(df_criterion) == expected_rows, \
    f"Expected {expected_rows} rows (14 tools × 28 criteria), got {len(df_criterion)}"
print(f"  ✓ Exactly {expected_rows} rows (14 tools × 28 criteria)")

# Validation: Check for duplicates
duplicates = df_criterion.duplicated(subset=['Tool', 'Criterion']).sum()
assert duplicates == 0, f"Found {duplicates} duplicate Tool × Criterion pairs"
print(f"  ✓ No duplicate Tool × Criterion pairs")

# Validation: Check for missing Jaccard values
missing_jaccard = df_criterion['Jaccard_Similarity'].isna().sum()
assert missing_jaccard == 0, f"Found {missing_jaccard} missing Jaccard values"
print(f"  ✓ No missing Jaccard values")

# Validation: Check Jaccard range
jaccard_values = df_criterion['Jaccard_Similarity']
assert (jaccard_values >= 0).all() and (jaccard_values <= 1).all(), \
    "Jaccard values outside [0, 1] range"
print(f"  ✓ All Jaccard values in [0, 1] range")

# Count tools and criteria
n_tools = df_criterion['Tool'].nunique()
n_criteria = df_criterion['Criterion'].nunique()
print(f"  Tools: {n_tools}")
print(f"  Criteria: {n_criteria}")

assert n_tools == 14, f"Expected 14 tools, found {n_tools}"
assert n_criteria == 28, f"Expected 28 criteria, found {n_criteria}"

# Validation: Every tool has exactly 28 criteria
tool_criterion_counts = df_criterion.groupby('Tool').size()
assert (tool_criterion_counts == 28).all(), \
    "Not all tools have exactly 28 criteria"
print(f"  ✓ Every tool has exactly 28 criteria")

# Pivot to wide format
print("\n3. Pivoting criterion Jaccard to wide format...")
df_wide = df_criterion.pivot(index='Tool', columns='Criterion', values='Jaccard_Similarity')

# Reorder columns to match CRITERIA_ORDER
df_wide = df_wide[CRITERIA_ORDER]
print(f"  ✓ Pivoted to wide format: {df_wide.shape}")
print(f"  ✓ Columns ordered: {len(CRITERIA_ORDER)} criteria")

# Load mean criterion Jaccard
print("\n4. Loading raw mean criterion Jaccard...")
df_mean = pd.read_csv(source_files['mean_jaccard'])
print(f"  Loaded: {len(df_mean)} rows")

# Use only raw_mean_jaccard
if 'raw_mean_jaccard' not in df_mean.columns:
    raise ValueError("Column 'raw_mean_jaccard' not found")

df_mean = df_mean[['tool', 'raw_mean_jaccard']].copy()
df_mean.rename(columns={
    'tool': 'Tool',
    'raw_mean_jaccard': 'Mean_Criterion_Jaccard'
}, inplace=True)

# Validation: Check range
assert (df_mean['Mean_Criterion_Jaccard'] >= 0).all() and \
       (df_mean['Mean_Criterion_Jaccard'] <= 1).all(), \
    "Mean Jaccard values outside [0, 1] range"
print(f"  ✓ Using raw_mean_jaccard")
print(f"  ✓ All values in [0, 1] range")

# Load total criterion Jaccard
print("\n5. Loading raw total criterion Jaccard...")
df_total = pd.read_csv(source_files['total_jaccard'])
print(f"  Loaded: {len(df_total)} rows")

# Use only raw_total_criterion_jaccard
if 'raw_total_criterion_jaccard' not in df_total.columns:
    raise ValueError("Column 'raw_total_criterion_jaccard' not found")

df_total = df_total[['tool', 'raw_total_criterion_jaccard']].copy()
df_total.rename(columns={
    'tool': 'Tool',
    'raw_total_criterion_jaccard': 'Total_Criterion_Jaccard'
}, inplace=True)

# Validation: Check range
assert (df_total['Total_Criterion_Jaccard'] >= 0).all() and \
       (df_total['Total_Criterion_Jaccard'] <= 1).all(), \
    "Total Jaccard values outside [0, 1] range"
print(f"  ✓ Using raw_total_criterion_jaccard")
print(f"  ✓ All values in [0, 1] range")

# Load 3-class Accuracy
print("\n6. Loading full-cohort 3-class Accuracy...")
df_accuracy = pd.read_csv(source_files['accuracy'])
print(f"  Loaded: {len(df_accuracy)} rows")

# Validation: Ensure we're using accuracy_3class, NOT Weighted F1
if 'Weighted F1 (3-class)' in df_accuracy.columns:
    raise ValueError("Weighted F1 column found - should not be used")
if 'accuracy_3class' not in df_accuracy.columns:
    raise ValueError("Column 'accuracy_3class' not found")

df_accuracy = df_accuracy[['tool', 'accuracy_3class']].copy()
df_accuracy.rename(columns={
    'tool': 'Tool',
    'accuracy_3class': 'Accuracy_3class'
}, inplace=True)

# Validation: Check range
assert (df_accuracy['Accuracy_3class'] >= 0).all() and \
       (df_accuracy['Accuracy_3class'] <= 1).all(), \
    "Accuracy values outside [0, 1] range"
print(f"  ✓ Using accuracy_3class (NOT Weighted F1)")
print(f"  ✓ All values in [0, 1] range")

# Merge all sources
print("\n7. Merging all sources by Tool...")

# Reset index to make Tool a column
df_wide = df_wide.reset_index()

# Merge criterion Jaccard with mean Jaccard
df_master = df_wide.merge(df_mean, on='Tool', how='inner')
print(f"  After merging mean Jaccard: {df_master.shape}")

# Merge with total Jaccard
df_master = df_master.merge(df_total, on='Tool', how='inner')
print(f"  After merging total Jaccard: {df_master.shape}")

# Merge with Accuracy
df_master = df_master.merge(df_accuracy, on='Tool', how='inner')
print(f"  After merging Accuracy: {df_master.shape}")

# Final validation
print("\n8. Final validation...")
print("-" * 80)

# Validation 1: Exactly 14 rows
assert len(df_master) == 14, f"Expected 14 rows, got {len(df_master)}"
print(f"  ✓ Exactly 14 rows (tool instances)")

# Validation 2: Exactly 32 columns (1 Tool + 28 criteria + 3 aggregates)
expected_cols = 1 + 28 + 3
assert len(df_master.columns) == expected_cols, \
    f"Expected {expected_cols} columns, got {len(df_master.columns)}"
print(f"  ✓ Exactly {expected_cols} columns (1 Tool + 28 criteria + 3 aggregates)")

# Validation 3: No missing values
missing_total = df_master.isna().sum().sum()
assert missing_total == 0, f"Found {missing_total} missing values"
print(f"  ✓ No missing values")

# Validation 4: Tool names match across sources
tools_criterion = set(df_wide['Tool'])
tools_mean = set(df_mean['Tool'])
tools_total = set(df_total['Tool'])
tools_accuracy = set(df_accuracy['Tool'])

assert tools_criterion == tools_mean == tools_total == tools_accuracy, \
    "Tool names do not match across sources"
print(f"  ✓ Tool names match across all sources")

# Validation 5: Verify column order
expected_columns = ['Tool'] + CRITERIA_ORDER + \
                  ['Mean_Criterion_Jaccard', 'Total_Criterion_Jaccard', 'Accuracy_3class']
assert list(df_master.columns) == expected_columns, \
    "Column order does not match expected order"
print(f"  ✓ Columns in correct order")

# Save master table
output_path = output_dir / 'correlation_master_table.csv'
df_master.to_csv(output_path, index=False)
print(f"\n✓ Saved: {output_path}")

# Print summary
print("\n" + "=" * 80)
print("CORRELATION MASTER TABLE SUMMARY")
print("=" * 80)

print(f"\nFinal Dimensions: {df_master.shape[0]} rows × {df_master.shape[1]} columns")

print(f"\nColumn Names ({len(df_master.columns)} total):")
print(f"  1. Tool")
print(f"  2-29. Criteria (28): {', '.join(CRITERIA_ORDER)}")
print(f"  30. Mean_Criterion_Jaccard")
print(f"  31. Total_Criterion_Jaccard")
print(f"  32. Accuracy_3class")

print(f"\nAll 14 Tool Names:")
for i, tool in enumerate(sorted(df_master['Tool']), 1):
    print(f"  {i:2d}. {tool}")

print(f"\nAggregate Metrics for All Tools:")
print("-" * 80)
df_agg = df_master[['Tool', 'Mean_Criterion_Jaccard', 
                    'Total_Criterion_Jaccard', 'Accuracy_3class']].copy()
df_agg = df_agg.sort_values('Accuracy_3class', ascending=False)

print(f"{'Tool':<20} {'Mean_Jaccard':>13} {'Total_Jaccard':>14} {'Accuracy':>10}")
print("-" * 80)
for _, row in df_agg.iterrows():
    print(f"{row['Tool']:<20} {row['Mean_Criterion_Jaccard']:>13.4f} "
          f"{row['Total_Criterion_Jaccard']:>14.4f} {row['Accuracy_3class']:>10.4f}")

print(f"\nCriterion-Specific Jaccard Statistics:")
print("-" * 80)
print(f"{'Criterion':<8} {'Min':>8} {'Max':>8} {'N_Zeros':>10} {'Mean':>8}")
print("-" * 80)

for criterion in CRITERIA_ORDER:
    values = df_master[criterion]
    n_zeros = (values == 0).sum()
    print(f"{criterion:<8} {values.min():>8.4f} {values.max():>8.4f} "
          f"{n_zeros:>10d} {values.mean():>8.4f}")

print(f"\nValidation Results:")
print("-" * 80)
print(f"  ✓ Exactly 14 rows / tool instances")
print(f"  ✓ Exactly 28 criterion columns")
print(f"  ✓ Criterion source contains exactly 392 rows = 14 × 28")
print(f"  ✓ Every tool has exactly one value for every criterion")
print(f"  ✓ No duplicated Tool × Criterion pairs")
print(f"  ✓ No missing criterion Jaccard values")
print(f"  ✓ All criterion Jaccards within [0, 1]")
print(f"  ✓ Mean, Total Jaccard and Accuracy all within [0, 1]")
print(f"  ✓ No missing values in final matrix")
print(f"  ✓ Tool names match across all sources")
print(f"  ✓ Accuracy from accuracy_3class (NOT Weighted F1)")

print(f"\nSource Files Used:")
print("-" * 80)
for name, path in source_files.items():
    print(f"  {name}:")
    print(f"    {path}")

print("\n" + "=" * 80)
print("MASTER TABLE BUILD COMPLETE")
print("=" * 80)
print("\nReady for correlation analysis.")
print("Do NOT calculate correlations or generate figures yet.")
print("=" * 80)
