#!/usr/bin/env python3
"""
Calculate Pairwise Pearson Correlations

GLOBAL EXPLORATORY ANALYSIS: Computes complete pairwise Pearson product-moment
correlations among 31 numeric variables (28 criteria + 3 aggregate metrics)
using 14 tool instances as observations.

This is parallel to the Spearman analysis but uses Pearson r instead of rho.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests

print("=" * 80)
print("CALCULATING PAIRWISE PEARSON CORRELATIONS")
print("=" * 80)

# Define paths
input_dir = Path(__file__).parent
output_dir = input_dir
master_table_path = input_dir / 'correlation_master_table.csv'

# Verify master table exists
if not master_table_path.exists():
    raise FileNotFoundError(f"Master table not found: {master_table_path}")

print(f"\n1. Loading master table...")
print(f"   {master_table_path}")

# Load master table
df = pd.read_csv(master_table_path)
print(f"   Loaded: {df.shape[0]} rows × {df.shape[1]} columns")

# Extract numeric variables (exclude Tool column)
numeric_cols = [col for col in df.columns if col != 'Tool']
df_numeric = df[numeric_cols].copy()

n_obs = len(df_numeric)
n_vars = len(numeric_cols)

print(f"\n2. Data summary...")
print(f"   Observations (tools): {n_obs}")
print(f"   Numeric variables: {n_vars}")
print(f"   Expected: 28 criteria + 3 aggregates = 31")

assert n_vars == 31, f"Expected 31 numeric variables, found {n_vars}"
assert n_obs == 14, f"Expected 14 observations, found {n_obs}"

# Identify ACMG/AMP criteria
criteria = [col for col in numeric_cols if col not in 
           ['Mean_Criterion_Jaccard', 'Total_Criterion_Jaccard', 'Accuracy_3class']]
aggregates = ['Mean_Criterion_Jaccard', 'Total_Criterion_Jaccard', 'Accuracy_3class']

print(f"   Criteria: {len(criteria)}")
print(f"   Aggregates: {len(aggregates)}")

assert len(criteria) == 28, f"Expected 28 criteria, found {len(criteria)}"
assert len(aggregates) == 3, f"Expected 3 aggregates, found {len(aggregates)}"

# Detect zero-variance variables
print(f"\n3. Detecting zero-variance variables...")
zero_var_vars = []
for col in numeric_cols:
    if df_numeric[col].std() == 0:
        zero_var_vars.append(col)
        print(f"   ✓ Zero variance: {col} (all values = {df_numeric[col].iloc[0]:.4f})")

if not zero_var_vars:
    print(f"   No zero-variance variables detected")
else:
    print(f"   Total zero-variance variables: {len(zero_var_vars)}")

# Verify expected zero-variance criteria
expected_zero_var = {'PP5', 'BP3', 'BP6'}
actual_zero_var = set(zero_var_vars)
assert actual_zero_var == expected_zero_var, \
    f"Expected zero-variance: {expected_zero_var}, got {actual_zero_var}"

# Calculate number of possible pairs
n_total_pairs = n_vars * (n_vars - 1) // 2  # Unique pairs (excluding diagonal)
print(f"\n4. Pairwise correlation setup...")
print(f"   Total unique pairs possible: {n_total_pairs}")

# Calculate Pearson correlations
print(f"\n5. Calculating Pearson product-moment correlations...")

# Initialize matrices
r_matrix = pd.DataFrame(np.nan, index=numeric_cols, columns=numeric_cols)
pval_matrix = pd.DataFrame(np.nan, index=numeric_cols, columns=numeric_cols)

# Store pairwise results for long format
pairwise_results = []

# Calculate correlations
n_valid_tests = 0
n_zero_var_excluded = 0

for i, var1 in enumerate(numeric_cols):
    for j, var2 in enumerate(numeric_cols):
        if i == j:
            # Diagonal: self-correlation
            if var1 not in zero_var_vars:
                r_matrix.loc[var1, var2] = 1.0
                pval_matrix.loc[var1, var2] = 0.0
            # else: leave as NA for zero-variance variables
        elif i < j:
            # Upper triangle: calculate correlation
            if var1 in zero_var_vars or var2 in zero_var_vars:
                # Zero-variance variable: correlation is undefined
                # Leave as NA, do not include in FDR correction
                r_matrix.loc[var1, var2] = np.nan
                r_matrix.loc[var2, var1] = np.nan
                pval_matrix.loc[var1, var2] = np.nan
                pval_matrix.loc[var2, var1] = np.nan
                
                # Store in pairwise results with NA
                pairwise_results.append({
                    'Variable_1': var1,
                    'Variable_2': var2,
                    'Pearson_r': np.nan,
                    'P_value': np.nan,
                    'n': n_obs
                })
                n_zero_var_excluded += 1
            else:
                # Valid correlation test
                r, pval = stats.pearsonr(df_numeric[var1], df_numeric[var2])
                
                # Store in matrices (symmetric)
                r_matrix.loc[var1, var2] = r
                r_matrix.loc[var2, var1] = r
                pval_matrix.loc[var1, var2] = pval
                pval_matrix.loc[var2, var1] = pval
                
                # Store in pairwise results
                pairwise_results.append({
                    'Variable_1': var1,
                    'Variable_2': var2,
                    'Pearson_r': r,
                    'P_value': pval,
                    'n': n_obs
                })
                n_valid_tests += 1

print(f"   Valid correlation tests: {n_valid_tests}")
print(f"   Excluded (zero variance): {n_zero_var_excluded}")
print(f"   Total pairs processed: {n_valid_tests + n_zero_var_excluded}")

# Apply Benjamini-Hochberg FDR correction
print(f"\n6. Applying Benjamini-Hochberg FDR correction...")

# Create DataFrame from pairwise results
df_pairwise = pd.DataFrame(pairwise_results)

# Separate valid tests from NA tests
valid_mask = ~df_pairwise['P_value'].isna()
df_valid = df_pairwise[valid_mask].copy()
df_na = df_pairwise[~valid_mask].copy()

print(f"   Valid tests for FDR correction: {len(df_valid)}")
print(f"   Excluded from FDR (zero variance): {len(df_na)}")

# Apply FDR correction only to valid tests
if len(df_valid) > 0:
    reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(
        df_valid['P_value'], 
        alpha=0.05, 
        method='fdr_bh'
    )
    df_valid['FDR_q_value'] = pvals_corrected
    df_valid['FDR_significant'] = reject
else:
    df_valid['FDR_q_value'] = []
    df_valid['FDR_significant'] = []

# For NA tests, set FDR q-value to NA
df_na['FDR_q_value'] = np.nan
df_na['FDR_significant'] = False

# Combine back
df_pairwise = pd.concat([df_valid, df_na], ignore_index=True)

# Create FDR matrix
fdr_matrix = pd.DataFrame(np.nan, index=numeric_cols, columns=numeric_cols)

# Fill FDR matrix
for _, row in df_pairwise.iterrows():
    var1 = row['Variable_1']
    var2 = row['Variable_2']
    fdr_q = row['FDR_q_value']
    
    fdr_matrix.loc[var1, var2] = fdr_q
    fdr_matrix.loc[var2, var1] = fdr_q

# Diagonal FDR values
for var in numeric_cols:
    if var not in zero_var_vars:
        fdr_matrix.loc[var, var] = 0.0

# Count significant results
n_raw_sig = (df_pairwise['P_value'] < 0.05).sum()
n_fdr_sig = (df_pairwise['FDR_significant'] == True).sum()

print(f"   Raw p < 0.05: {n_raw_sig}")
print(f"   FDR q < 0.05: {n_fdr_sig}")

# Validation
print(f"\n7. Validation...")
# Check r values in range
valid_r = df_pairwise['Pearson_r'].dropna()
assert (valid_r >= -1).all() and (valid_r <= 1).all(), \
    "Pearson r values outside [-1, 1] range"
print(f"   ✓ All valid Pearson r values in [-1, 1]")

valid_p = df_pairwise['P_value'].dropna()
assert (valid_p >= 0).all() and (valid_p <= 1).all(), \
    "p-values outside [0, 1] range"
print(f"   ✓ All valid p-values in [0, 1]")

valid_q = df_pairwise['FDR_q_value'].dropna()
assert (valid_q >= 0).all() and (valid_q <= 1).all(), \
    "FDR q-values outside [0, 1] range"
print(f"   ✓ All valid FDR q-values in [0, 1]")

# Save outputs
print(f"\n8. Saving outputs...")

# Save matrices
r_matrix_path = output_dir / 'pearson_r_matrix.csv'
pval_matrix_path = output_dir / 'pearson_pvalue_matrix.csv'
fdr_matrix_path = output_dir / 'pearson_fdr_matrix.csv'

r_matrix.to_csv(r_matrix_path)
pval_matrix.to_csv(pval_matrix_path)
fdr_matrix.to_csv(fdr_matrix_path)

print(f"   ✓ {r_matrix_path}")
print(f"   ✓ {pval_matrix_path}")
print(f"   ✓ {fdr_matrix_path}")

# Save pairwise long format
pairwise_path = output_dir / 'pearson_pairwise_long.csv'
df_pairwise.to_csv(pairwise_path, index=False)
print(f"   ✓ {pairwise_path}")

# Create QC summary
print(f"\n9. Creating QC summary...")

qc_summary = []
qc_summary.append("=" * 80)
qc_summary.append("PEARSON CORRELATION QC SUMMARY")
qc_summary.append("=" * 80)
qc_summary.append("")
qc_summary.append(f"Number of observations (tools): {n_obs}")
qc_summary.append(f"Number of numeric variables: {n_vars}")
qc_summary.append("")
qc_summary.append("Zero-variance variables:")
if zero_var_vars:
    for var in zero_var_vars:
        qc_summary.append(f"  - {var} (all values = {df_numeric[var].iloc[0]:.4f})")
else:
    qc_summary.append("  None detected")
qc_summary.append("")
qc_summary.append(f"Number of unique possible pairs: {n_total_pairs}")
qc_summary.append(f"Number of valid correlation tests: {n_valid_tests}")
qc_summary.append(f"Number excluded (zero variance): {n_zero_var_excluded}")
qc_summary.append("")
qc_summary.append(f"Number with raw p < 0.05: {n_raw_sig}")
qc_summary.append(f"Number with FDR q < 0.05: {n_fdr_sig}")
qc_summary.append("")
qc_summary.append("=" * 80)

qc_path = output_dir / 'pearson_correlation_qc_summary.txt'
with open(qc_path, 'w') as f:
    f.write('\n'.join(qc_summary))

print(f"   ✓ {qc_path}")

# Print aggregate metric correlations
print("\n" + "=" * 80)
print("KEY AGGREGATE METRIC CORRELATIONS")
print("=" * 80)

key_pairs = [
    ('Mean_Criterion_Jaccard', 'Accuracy_3class'),
    ('Total_Criterion_Jaccard', 'Accuracy_3class'),
    ('Mean_Criterion_Jaccard', 'Total_Criterion_Jaccard')
]

print(f"\n{'Variable_1':<25} {'Variable_2':<25} {'r':>8} {'p':>10} {'FDR_q':>10}")
print("-" * 80)

for var1, var2 in key_pairs:
    pair = df_pairwise[
        ((df_pairwise['Variable_1'] == var1) & (df_pairwise['Variable_2'] == var2)) |
        ((df_pairwise['Variable_1'] == var2) & (df_pairwise['Variable_2'] == var1))
    ]
    
    if len(pair) > 0:
        row = pair.iloc[0]
        r = row['Pearson_r']
        p = row['P_value']
        q = row['FDR_q_value']
        
        print(f"{var1:<25} {var2:<25} {r:>8.4f} {p:>10.4f} {q:>10.4f}")

# Final summary
print("\n" + "=" * 80)
print("GLOBAL PEARSON CORRELATION CALCULATION COMPLETE")
print("=" * 80)

print(f"\nOutputs created:")
print(f"  ✓ pearson_r_matrix.csv (31 × 31)")
print(f"  ✓ pearson_pvalue_matrix.csv (31 × 31)")
print(f"  ✓ pearson_fdr_matrix.csv (31 × 31)")
print(f"  ✓ pearson_pairwise_long.csv ({len(df_pairwise)} pairs)")
print(f"  ✓ pearson_correlation_qc_summary.txt")

print(f"\nThis is a GLOBAL exploratory analysis.")
print(f"For PRIMARY criterion-aggregate analysis, run:")
print(f"  calculate_primary_criterion_macro_pearson_correlations.py")

print("\n" + "=" * 80)
