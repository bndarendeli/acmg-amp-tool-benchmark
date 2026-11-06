#!/usr/bin/env python3
"""
Calculate Pairwise Pearson Correlations (Figure 5)

Observations are the 14 ClinGen tool instances in `correlation_master_table.csv`
(1 Tool column + 28 criterion-level Jaccard variables + 3 aggregate metrics).
The master table is read only; it is never modified or filtered by row.

Analysis design
---------------
1. Criterion sparsity QC (`criterion_correlation_qc.csv`).
   A zero criterion-level Jaccard means zero overlap between the tool's and the
   reference's variant sets for that criterion; it does NOT imply that the tool
   never assigned the criterion.
2. Primary eligibility rule (exploratory, not a prespecified or established cutoff):
   strictly positive Jaccard in >= 5 of the 14 tool instances and nonzero variance.
   All 14 observations (including zeros) are kept for every eligible criterion.
   The 3/14, 4/14 and 7/14 thresholds are evaluated in sensitivity analyses.
3. Two Benjamini-Hochberg FDR families:
   Family A - all unique criterion-criterion pairs among eligible criteria.
   Family B - eligible criterion x aggregate pairs plus the three
              aggregate-aggregate pairs (MJ, TJ, 3-class accuracy).
   Aggregate metrics are used exactly as supplied (not recomputed after filtering).
4. Sensitivity analyses (written to `sensitivity/`):
   minimum positive-Jaccard thresholds (3, 4, 5, 7 of 14), leave-one-tool-instance-out
   and one-representative-per-related-tool-family subsets (for the 5/14 primary and the
   7/14 comparison eligibility sets; files tagged `_min5of14` / `_min7of14`), and a reference reproduction
   of the zero-variance-only (378-test) sensitivity analysis.

Pearson p-values assume independent, bivariate-normal observations; with 14
partly related tool instances and zero-inflated variables they are descriptive.
"""

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

print("=" * 80)
print("CALCULATING PAIRWISE PEARSON CORRELATIONS (FIGURE 5)")
print("=" * 80)

input_dir = Path(__file__).resolve().parent
output_dir = input_dir
sens_dir = output_dir / 'sensitivity'
sens_dir.mkdir(parents=True, exist_ok=True)
base_dir = input_dir.parent
master_table_path = input_dir / 'correlation_master_table.csv'

SOURCE_FILES = {
    'criterion_jaccard': base_dir / 'fig3' / 'jaccard_per_criterion_clingen_28012026.csv',
    'mean_jaccard': base_dir / 'fig6' / 'outputs' / 'fig6_capability_adjusted_summary.csv',
    'total_jaccard': base_dir / 'fig6' / 'outputs' / 'clingen_total_criterion_jaccard.csv',
    'accuracy': base_dir / 'fig4' / 'data' / 'reconstructed' / 'ClinGen_reconstructed_metrics.csv',
}

CRITERIA_ORDER = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]
AGGREGATES = ['Mean_Criterion_Jaccard', 'Total_Criterion_Jaccard', 'Accuracy_3class']
AGGREGATE_LABELS = {'Mean_Criterion_Jaccard': 'MJ', 'Total_Criterion_Jaccard': 'TJ',
                    'Accuracy_3class': '3CA'}

PRIMARY_MIN_POSITIVE = 5   # exploratory threshold; 3/14 and 7/14 evaluated in sensitivity analyses
SENSITIVITY_THRESHOLDS = [3, 4, 5, 7]
ROBUSTNESS_THRESHOLDS = [PRIMARY_MIN_POSITIVE, 7]   # leave-one-out / tool-family analyses
ALPHA = 0.05
EXPECTED_PRIMARY_ELIGIBLE = ['PVS1', 'PS1', 'PS3', 'PM1', 'PM2', 'PM4', 'PM5', 'PP2', 'PP3',
                             'BA1', 'BS1', 'BS2', 'BP1', 'BP4', 'BP7']
EXPECTED_FAMILY_SIZES = {'A': 105, 'B': 48}

FAMILY_A = 'A_criterion_criterion'
FAMILY_B = 'B_criterion_aggregate_and_aggregate_aggregate'

# Related tool instances (shared code base / versions / deployments).
# Primary definition: the two documented same-tool pairs.
# Extended definition: additionally treats AutoGVP as InterVar-derived (AutoGVP
# runs a modified InterVar; Bioinformatics 2024, doi:10.1093/bioinformatics/btae114).
TOOL_FAMILY_DEFINITIONS = {
    'related_tool_families': {
        'CharGer': ['CharGer_Local', 'CharGer_Online'],
        'InterVar': ['InterVar_2018', 'InterVar_2025'],
    },
    'related_tool_families_extended': {
        'CharGer': ['CharGer_Local', 'CharGer_Online'],
        'InterVar_derived': ['InterVar_2018', 'InterVar_2025', 'AutoGVP'],
    },
}
# Pre-specified, data-independent representative rule:
#   most recent release of the parent tool where versions differ (InterVar_2025);
#   the locally executed standalone software where only the deployment differs
#   (CharGer_Local). Results for every possible representative combination are
#   also reported so the choice can be audited.
PRESPECIFIED_REPRESENTATIVE = {'CharGer': 'CharGer_Local', 'InterVar': 'InterVar_2025',
                               'InterVar_derived': 'InterVar_2025'}
FOCUS_INSTANCES = ['Genebe', 'Franklin']


def label(var):
    return AGGREGATE_LABELS.get(var, var)


def pair_type(v1, v2):
    n_agg = (v1 in AGGREGATES) + (v2 in AGGREGATES)
    return ['criterion-criterion', 'criterion-aggregate', 'aggregate-aggregate'][n_agg]


def pearson(x, y):
    """Pearson r and two-sided p; NaN when either variable has zero variance."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, np.nan
    r, p = stats.pearsonr(x, y)
    return float(r), float(p)


def apply_bh(df, family_col='FDR_family'):
    """BH-FDR within each family, over valid (non-NaN) tests only."""
    df = df.copy()
    df['FDR_q_value'] = np.nan
    df['FDR_family_size'] = 0
    for fam, idx in df.groupby(family_col).groups.items():
        sub = df.loc[idx]
        valid = sub.index[sub['P_value'].notna()]
        if len(valid):
            df.loc[valid, 'FDR_q_value'] = multipletests(df.loc[valid, 'P_value'], method='fdr_bh')[1]
        df.loc[idx, 'FDR_family_size'] = len(valid)
    df['FDR_significant'] = df['FDR_q_value'].lt(ALPHA)  # NaN -> False
    return df


def correlate(df, criteria, include_family_b=True, aggregates=AGGREGATES):
    """Family A (criterion pairs) and optionally Family B; BH within family."""
    rows = []
    for v1, v2 in itertools.combinations(criteria, 2):
        rows.append((v1, v2, FAMILY_A))
    if include_family_b:
        for c in criteria:
            for a in aggregates:
                rows.append((c, a, FAMILY_B))
        for a1, a2 in itertools.combinations(aggregates, 2):
            rows.append((a1, a2, FAMILY_B))
    out = []
    for v1, v2, fam in rows:
        r, p = pearson(df[v1], df[v2])
        out.append({'Variable_1': v1, 'Variable_2': v2, 'Pair_type': pair_type(v1, v2),
                    'FDR_family': fam, 'n': len(df), 'Pearson_r': r, 'P_value': p})
    return apply_bh(pd.DataFrame(out))


def positive_counts(df, criteria):
    return {c: int((df[c] > 0).sum()) for c in criteria}


def eligible_criteria(df, criteria, min_positive):
    return [c for c in criteria if (df[c] > 0).sum() >= min_positive and df[c].std() > 0]


def find_pair(long_df, v1, v2):
    m = (((long_df['Variable_1'] == v1) & (long_df['Variable_2'] == v2)) |
         ((long_df['Variable_1'] == v2) & (long_df['Variable_2'] == v1)))
    sub = long_df[m]
    return sub.iloc[0] if len(sub) else None


def extreme_pair(long_df, positive=True):
    sub = long_df.dropna(subset=['Pearson_r'])
    if positive:
        sub = sub[sub['Pearson_r'] > 0]
        row = sub.loc[sub['Pearson_r'].idxmax()] if len(sub) else None
    else:
        sub = sub[sub['Pearson_r'] < 0]
        row = sub.loc[sub['Pearson_r'].idxmin()] if len(sub) else None
    if row is None:
        return 'none', np.nan, np.nan
    return f"{row['Variable_1']}-{row['Variable_2']}", row['Pearson_r'], row['FDR_q_value']


# ---------------------------------------------------------------------------
# 1. Load and validate the master table (read only)
# ---------------------------------------------------------------------------
print("\n1. Loading and validating master table (read only)...")
if not master_table_path.exists():
    raise FileNotFoundError(f"Master table not found: {master_table_path}")
df = pd.read_csv(master_table_path)
assert df.shape == (14, 32), f"Expected 14 rows x 32 columns, got {df.shape}"
assert list(df.columns) == ['Tool'] + CRITERIA_ORDER + AGGREGATES, "Unexpected column order"
assert df['Tool'].is_unique and not df.isna().any().any()
tools = list(df['Tool'])
n_obs = len(df)
print(f"   ✓ {df.shape[0]} rows x {df.shape[1]} columns; tool instances: {', '.join(tools)}")

# Source values unchanged: master table must equal the Fig3/Fig4/Fig6 sources
src = pd.read_csv(SOURCE_FILES['criterion_jaccard']).pivot(
    index='Tool', columns='Criterion', values='Jaccard_Similarity')[CRITERIA_ORDER]
src = src.loc[tools]
assert np.allclose(src.values, df[CRITERIA_ORDER].values, rtol=0, atol=1e-12), \
    "Criterion Jaccard values differ from Figure 3 source"
agg_src = {
    'Mean_Criterion_Jaccard': pd.read_csv(SOURCE_FILES['mean_jaccard']).set_index('tool')['raw_mean_jaccard'],
    'Total_Criterion_Jaccard': pd.read_csv(SOURCE_FILES['total_jaccard']).set_index('tool')['raw_total_criterion_jaccard'],
    'Accuracy_3class': pd.read_csv(SOURCE_FILES['accuracy']).set_index('tool')['accuracy_3class'],
}
for col, s in agg_src.items():
    assert np.allclose(s.loc[tools].values, df[col].values, rtol=0, atol=1e-12), \
        f"{col} differs from its source file"
print("   ✓ Criterion Jaccard and aggregate values identical to source files")
data = df.set_index('Tool')

# ---------------------------------------------------------------------------
# 2. Criterion sparsity QC
# ---------------------------------------------------------------------------
print("\n2. Criterion sparsity QC...")
qc_rows = []
for c in CRITERIA_ORDER:
    v = data[c]
    n_pos = int((v > 0).sum())
    sd = float(v.std(ddof=1))
    reasons = []
    if sd == 0:
        reasons.append(f"zero variance (all {n_obs} values = {v.iloc[0]:g})")
    if n_pos < PRIMARY_MIN_POSITIVE:
        reasons.append(f"positive Jaccard in {n_pos}/{n_obs} instances "
                       f"(< {PRIMARY_MIN_POSITIVE}/{n_obs} required)")
    row = {
        'Criterion': c,
        'N_tool_instances': int(v.notna().sum()),
        'N_zero_jaccard': int((v == 0).sum()),
        'N_positive_jaccard': n_pos,
        'Proportion_positive': n_pos / n_obs,
        'N_distinct_values': int(v.nunique()),
        'SD_jaccard': sd,
        'Mean_jaccard': float(v.mean()),
        'Max_jaccard': float(v.max()),
        'Instances_with_positive_jaccard': ';'.join(v.index[v > 0]),
        'Eligible_primary': not reasons,
        'Exclusion_reason': '; '.join(reasons),
    }
    for t in SENSITIVITY_THRESHOLDS:
        row[f'Eligible_min{t}of{n_obs}'] = bool(n_pos >= t and sd > 0)
    qc_rows.append(row)
qc = pd.DataFrame(qc_rows)
qc.to_csv(output_dir / 'criterion_correlation_qc.csv', index=False)

eligible = list(qc.loc[qc['Eligible_primary'], 'Criterion'])
excluded = list(qc.loc[~qc['Eligible_primary'], 'Criterion'])
assert eligible == eligible_criteria(data, CRITERIA_ORDER, PRIMARY_MIN_POSITIVE)
assert eligible == EXPECTED_PRIMARY_ELIGIBLE, \
    f"Eligible criteria differ from expectation: {eligible}"
for c in ['PS2', 'PM3', 'PP1', 'PP4']:
    assert c in excluded, f"{c} must be excluded from the primary analysis"
print(f"   Eligible ({len(eligible)}): {', '.join(eligible)}")
print(f"   Excluded ({len(excluded)}): {', '.join(excluded)}")
print("   ✓ Eligible list verified programmatically")

# ---------------------------------------------------------------------------
# 3. Primary analysis: Family A and Family B
# ---------------------------------------------------------------------------
print("\n3. Primary correlation analysis...")
primary = correlate(data, eligible)
primary.insert(0, 'Analysis', f'primary_min{PRIMARY_MIN_POSITIVE}of{n_obs}')
n_a = int((primary['FDR_family'] == FAMILY_A).sum())
n_b = int((primary['FDR_family'] == FAMILY_B).sum())
assert n_a == len(eligible) * (len(eligible) - 1) // 2 == EXPECTED_FAMILY_SIZES['A'], f"Family A has {n_a} tests"
assert n_b == len(eligible) * len(AGGREGATES) + 3 == EXPECTED_FAMILY_SIZES['B'], f"Family B has {n_b} tests"
assert primary['Pearson_r'].notna().all(), "Undefined correlation in primary analysis"
assert not primary.loc[primary['Pearson_r'].isna(), 'FDR_significant'].any()
print(f"   Family A ({FAMILY_A}): {n_a} tests, "
      f"{int(primary.loc[primary.FDR_family == FAMILY_A, 'FDR_significant'].sum())} q < {ALPHA}")
print(f"   Family B ({FAMILY_B}): {n_b} tests, "
      f"{int(primary.loc[primary.FDR_family == FAMILY_B, 'FDR_significant'].sum())} q < {ALPHA}")

variables = eligible + AGGREGATES
r_matrix = pd.DataFrame(np.nan, index=variables, columns=variables)
p_matrix = r_matrix.copy()
q_matrix = r_matrix.copy()
fam_matrix = pd.DataFrame('', index=variables, columns=variables)
for _, row in primary.iterrows():
    for a, b in [(row.Variable_1, row.Variable_2), (row.Variable_2, row.Variable_1)]:
        r_matrix.loc[a, b] = row.Pearson_r
        p_matrix.loc[a, b] = row.P_value
        q_matrix.loc[a, b] = row.FDR_q_value
        fam_matrix.loc[a, b] = row.FDR_family
for v in variables:
    r_matrix.loc[v, v] = 1.0
    fam_matrix.loc[v, v] = 'diagonal'
# Diagonal p and q are left empty (self-correlations are not tests)

for m in [r_matrix, p_matrix, q_matrix]:
    vals = m.values.astype(float)
    off = ~np.eye(len(variables), dtype=bool)
    assert np.allclose(vals[off], vals.T[off], rtol=0, atol=0), "Matrix not symmetric"
    assert not np.isnan(vals[off]).any()
assert (fam_matrix.values == fam_matrix.values.T).all()
assert not set(excluded) & set(r_matrix.index), "Excluded criteria present in primary matrix"
print(f"   ✓ Matrices {len(variables)}x{len(variables)} symmetric; excluded criteria absent")

r_matrix.to_csv(output_dir / 'pearson_r_matrix.csv')
p_matrix.to_csv(output_dir / 'pearson_pvalue_matrix.csv')
q_matrix.to_csv(output_dir / 'pearson_fdr_matrix.csv')
fam_matrix.to_csv(output_dir / 'pearson_fdr_family_matrix.csv')

# ---------------------------------------------------------------------------
# 4. Sensitivity: zero-variance-only analysis (378 tests, one family)
# ---------------------------------------------------------------------------
print("\n4. Sensitivity analysis (zero-variance exclusion only)...")
nonzero_var = [c for c in CRITERIA_ORDER if data[c].std() > 0]
nonzero_variance_sensitivity_vars = nonzero_var + AGGREGATES
nonzero_variance_sensitivity = []
for v1, v2 in itertools.combinations(nonzero_variance_sensitivity_vars, 2):
    r, p = pearson(data[v1], data[v2])
    nonzero_variance_sensitivity.append({'Variable_1': v1, 'Variable_2': v2, 'Pair_type': pair_type(v1, v2),
                   'FDR_family': 'nonzero_variance_sensitivity_single_family_all_pairs', 'n': n_obs,
                   'Pearson_r': r, 'P_value': p})
nonzero_variance_sensitivity = apply_bh(pd.DataFrame(nonzero_variance_sensitivity))
nonzero_variance_sensitivity.insert(0, 'Analysis', 'nonzero_variance_sensitivity_zero_variance_exclusion_only')
assert len(nonzero_variance_sensitivity) == 378
nonzero_variance_sensitivity.to_csv(sens_dir / 'sensitivity_nonzero_variance_pairwise_long.csv', index=False)
print(f"   378 tests; q < {ALPHA}: {int(nonzero_variance_sensitivity['FDR_significant'].sum())}")

# ---------------------------------------------------------------------------
# 5. Sensitivity: minimum positive-Jaccard thresholds
# ---------------------------------------------------------------------------
print("\n5. Threshold sensitivity analysis...")
nonzero_variance_sensitivity_cc = nonzero_variance_sensitivity[(nonzero_variance_sensitivity.Pair_type == 'criterion-criterion') & nonzero_variance_sensitivity.FDR_significant]
nonzero_variance_sensitivity_top = nonzero_variance_sensitivity_cc.reindex(nonzero_variance_sensitivity_cc['Pearson_r'].abs().sort_values(ascending=False).index).head(10)
def top_sig(res, k=5):
    a = res[(res.FDR_family == FAMILY_A) & res.FDR_significant]
    return a.reindex(a['Pearson_r'].abs().sort_values(ascending=False).index).head(k)


major_sources = [('nonzero_variance_sensitivity_top10_abs_r_FDR_sig', nonzero_variance_sensitivity_top)]
for t in sorted(set(ROBUSTNESS_THRESHOLDS), reverse=True):
    res_t = primary if t == PRIMARY_MIN_POSITIVE else correlate(data, eligible_criteria(data, CRITERIA_ORDER, t))
    major_sources.append((f'min{t}of{n_obs}{"_primary" if t == PRIMARY_MIN_POSITIVE else ""}'
                          '_top5_abs_r_FDR_sig', top_sig(res_t)))
major_pairs = []
for src_name, sub in major_sources:
    for _, row in sub.iterrows():
        key = (row.Variable_1, row.Variable_2)
        if key not in [(m[0], m[1]) for m in major_pairs]:
            major_pairs.append((row.Variable_1, row.Variable_2, src_name))

thr_summary, thr_long, thr_major = [], [], []
for t in SENSITIVITY_THRESHOLDS:
    crit_t = eligible_criteria(data, CRITERIA_ORDER, t)
    res = correlate(data, crit_t)
    res.insert(0, 'Analysis', f'threshold_min{t}of{n_obs}')
    thr_long.append(res)
    a = res[res.FDR_family == FAMILY_A]
    b = res[res.FDR_family == FAMILY_B]
    pos_pair, pos_r, pos_q = extreme_pair(a, True)
    neg_pair, neg_r, neg_q = extreme_pair(a, False)
    thr_summary.append({
        'Min_positive_jaccard': f'{t}/{n_obs}',
        'Is_primary': t == PRIMARY_MIN_POSITIVE,
        'N_eligible_criteria': len(crit_t),
        'Eligible_criteria': ';'.join(crit_t),
        'N_family_A_pairs_tested': int(a['P_value'].notna().sum()),
        'N_family_A_FDR_significant': int(a['FDR_significant'].sum()),
        'N_family_A_FDR_sig_positive': int((a['FDR_significant'] & (a['Pearson_r'] > 0)).sum()),
        'N_family_A_FDR_sig_negative': int((a['FDR_significant'] & (a['Pearson_r'] < 0)).sum()),
        'Highest_positive_pair': pos_pair, 'Highest_positive_r': pos_r, 'Highest_positive_q': pos_q,
        'Most_negative_pair': neg_pair, 'Most_negative_r': neg_r, 'Most_negative_q': neg_q,
        'N_family_B_tests': int(b['P_value'].notna().sum()),
        'N_family_B_FDR_significant': int(b['FDR_significant'].sum()),
    })
    for v1, v2, src_name in major_pairs:
        row = find_pair(a, v1, v2)
        thr_major.append({
            'Pair': f'{v1}-{v2}', 'Source_of_pair': src_name, 'Min_positive_jaccard': f'{t}/{n_obs}',
            'Both_eligible': row is not None,
            'Pearson_r': row.Pearson_r if row is not None else np.nan,
            'FDR_q_value': row.FDR_q_value if row is not None else np.nan,
            'FDR_significant': bool(row.FDR_significant) if row is not None else False,
        })
thr_summary = pd.DataFrame(thr_summary)
thr_summary.to_csv(sens_dir / 'sensitivity_threshold_summary.csv', index=False)
pd.concat(thr_long, ignore_index=True).to_csv(sens_dir / 'sensitivity_threshold_pairwise_long.csv', index=False)
pd.DataFrame(thr_major).to_csv(sens_dir / 'sensitivity_threshold_major_pairs.csv', index=False)
print(thr_summary[['Min_positive_jaccard', 'N_eligible_criteria', 'N_family_A_pairs_tested',
                   'N_family_A_FDR_significant']].to_string(index=False))

# ---------------------------------------------------------------------------
# 6-7. Robustness: leave-one-out and tool-instance dependence
#      Run for the primary threshold and for each comparison threshold;
#      outputs carry an explicit `_min{k}of14` tag.
# ---------------------------------------------------------------------------
print("\n6-7. Robustness analyses (leave-one-out; tool-instance dependence)...")
# Empirical similarity of criterion profiles between instances (descriptive)
sim_rows = []
for t1, t2 in itertools.combinations(tools, 2):
    x, y = data.loc[t1, CRITERIA_ORDER], data.loc[t2, CRITERIA_ORDER]
    sim_rows.append({
        'Instance_1': t1, 'Instance_2': t2,
        'Profile_pearson_r_28_criteria': pearson(x, y)[0],
        'N_identical_positive_criterion_jaccard': int(((x == y) & (x > 0)).sum()),
        'N_identical_criterion_jaccard_incl_zero': int((x == y).sum()),
    })
sim = pd.DataFrame(sim_rows).sort_values('Profile_pearson_r_28_criteria', ascending=False)

fam_rows = []
for def_name, fams in TOOL_FAMILY_DEFINITIONS.items():
    members = {m: f for f, ms in fams.items() for m in ms}
    for t in tools:
        fam_rows.append({'Family_definition': def_name, 'Tool_instance': t,
                         'Related_family': members.get(t, t),
                         'Prespecified_representative': members.get(t) is None or
                         PRESPECIFIED_REPRESENTATIVE[members[t]] == t})
fam_df = pd.DataFrame(fam_rows)
fam_df.to_csv(sens_dir / 'tool_instance_families.csv', index=False)
sim.to_csv(sens_dir / 'tool_instance_profile_similarity.csv', index=False)


def robustness(elig, base, tag):
    """Leave-one-out and one-representative-per-family analyses for `base` pairs."""
    # Leave-one-tool-instance-out (BH recomputed within each reduced family)
    loo_detail = []
    for dropped in tools:
        res = correlate(data.drop(index=dropped), elig)
        res.insert(0, 'Dropped_instance', dropped)
        loo_detail.append(res)
    loo_detail = pd.concat(loo_detail, ignore_index=True)
    loo_detail.to_csv(sens_dir / f'sensitivity_leave_one_out_detail_{tag}.csv', index=False)

    loo_rows = []
    for _, row in base.iterrows():
        d = loo_detail[(loo_detail.Variable_1 == row.Variable_1) & (loo_detail.Variable_2 == row.Variable_2)]
        rs = d.set_index('Dropped_instance')['Pearson_r']
        delta = (rs - row.Pearson_r).abs()
        out = {
            'Variable_1': row.Variable_1, 'Variable_2': row.Variable_2, 'Pair_type': row.Pair_type,
            'FDR_family': row.FDR_family, 'Full_r': row.Pearson_r, 'Full_p': row.P_value,
            'Full_q': row.FDR_q_value, 'Full_FDR_significant': row.FDR_significant,
            'LOO_min_r': rs.min(), 'LOO_max_r': rs.max(),
            'LOO_max_abs_change_r': delta.max(), 'Instance_with_max_change': delta.idxmax(),
            'N_LOO_undefined': int(rs.isna().sum()),
            'LOO_sign_always_preserved': bool((np.sign(rs) == np.sign(row.Pearson_r)).all()),
            'LOO_max_p': d['P_value'].max(),
            'LOO_all_nominal_p_lt_0.05': bool((d['P_value'] < ALPHA).all()),
            'LOO_all_FDR_q_lt_0.05': bool(d['FDR_significant'].all()),
        }
        for inst in FOCUS_INSTANCES:
            out[f'r_without_{inst}'] = rs.get(inst, np.nan)
            out[f'q_without_{inst}'] = d.set_index('Dropped_instance')['FDR_q_value'].get(inst, np.nan)
        both = data.drop(index=FOCUS_INSTANCES)
        out['r_without_' + '_and_'.join(FOCUS_INSTANCES)] = pearson(both[row.Variable_1], both[row.Variable_2])[0]
        loo_rows.append(out)
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(sens_dir / f'sensitivity_leave_one_out_summary_{tag}.csv', index=False)

    # One representative per related family: pre-specified and all combinations
    dep_long, dep_scen = [], []
    for def_name, fams in TOOL_FAMILY_DEFINITIONS.items():
        singletons = [t for t in tools if not any(t in ms for ms in fams.values())]
        fam_names = list(fams)
        for combo in itertools.product(*[fams[f] for f in fam_names]):
            retained = [t for t in tools if t in singletons or t in combo]
            is_pre = all(PRESPECIFIED_REPRESENTATIVE[f] == c for f, c in zip(fam_names, combo))
            scen = f"{def_name}:{'+'.join(combo)}"
            sub = data.loc[retained]
            res = correlate(sub, elig)
            res.insert(0, 'Scenario', scen)
            res.insert(1, 'Family_definition', def_name)
            res.insert(2, 'Prespecified', is_pre)
            dep_long.append(res)
            a = res[res.FDR_family == FAMILY_A]
            s = {'Scenario': scen, 'Family_definition': def_name, 'Prespecified': is_pre,
                 'Representatives': ';'.join(combo), 'N_instances': len(retained),
                 'Retained_instances': ';'.join(retained),
                 'Eligible_criteria_with_zero_variance': ';'.join(c for c in elig if sub[c].std() == 0),
                 'Min_positive_count_among_primary_criteria': min(positive_counts(sub, elig).values()),
                 'N_family_A_tests': int(a['P_value'].notna().sum()),
                 'N_family_A_FDR_significant': int(a['FDR_significant'].sum()),
                 'N_family_B_FDR_significant': int(res.loc[res.FDR_family == FAMILY_B, 'FDR_significant'].sum())}
            for v1, v2 in [('Mean_Criterion_Jaccard', 'Accuracy_3class'), ('Total_Criterion_Jaccard', 'Accuracy_3class')]:
                row = find_pair(res, v1, v2)
                k = f'{label(v1)}_vs_{label(v2)}'
                s[f'{k}_r'], s[f'{k}_p'], s[f'{k}_q'] = row.Pearson_r, row.P_value, row.FDR_q_value
            dep_scen.append(s)
    dep_long = pd.concat(dep_long, ignore_index=True)
    dep_scen = pd.DataFrame(dep_scen)
    dep_long.to_csv(sens_dir / f'sensitivity_tool_family_pairwise_long_{tag}.csv', index=False)
    dep_scen.to_csv(sens_dir / f'sensitivity_tool_family_scenarios_{tag}.csv', index=False)

    dep_pair_rows = []
    for _, row in base.iterrows():
        out = {'Variable_1': row.Variable_1, 'Variable_2': row.Variable_2, 'FDR_family': row.FDR_family,
               'Full_r': row.Pearson_r, 'Full_q': row.FDR_q_value}
        for def_name in TOOL_FAMILY_DEFINITIONS:
            d = dep_long[(dep_long.Family_definition == def_name) & (dep_long.Variable_1 == row.Variable_1)
                         & (dep_long.Variable_2 == row.Variable_2)]
            pre = d[d.Prespecified].iloc[0]
            out[f'{def_name}_prespecified_r'] = pre.Pearson_r
            out[f'{def_name}_prespecified_q'] = pre.FDR_q_value
            out[f'{def_name}_min_r'] = d['Pearson_r'].min()
            out[f'{def_name}_max_r'] = d['Pearson_r'].max()
            out[f'{def_name}_sign_preserved_all_combinations'] = bool((np.sign(d['Pearson_r']) == np.sign(row.Pearson_r)).all())
            out[f'{def_name}_FDR_sig_all_combinations'] = bool(d['FDR_significant'].all())
        dep_pair_rows.append(out)
    dep_pairs = pd.DataFrame(dep_pair_rows)
    dep_pairs.to_csv(sens_dir / f'sensitivity_tool_family_pair_summary_{tag}.csv', index=False)
    return loo, dep_scen, dep_pairs


def add_robust_flags(base, loo, dep_pairs):
    """Robust = BH q < 0.05 in its family AND same sign with nominal p < 0.05 in every
    leave-one-out subsample AND BH q < 0.05 (recomputed) in every
    one-representative-per-family combination (primary and extended definitions)."""
    merge_cols = ['Variable_1', 'Variable_2']
    out = base.merge(loo[merge_cols + ['LOO_min_r', 'LOO_max_r', 'LOO_max_abs_change_r',
                                       'Instance_with_max_change', 'LOO_sign_always_preserved',
                                       'LOO_all_nominal_p_lt_0.05', 'r_without_Genebe',
                                       'r_without_Franklin']], on=merge_cols)
    out = out.merge(dep_pairs[merge_cols + [
        'related_tool_families_prespecified_r', 'related_tool_families_prespecified_q',
        'related_tool_families_FDR_sig_all_combinations',
        'related_tool_families_extended_FDR_sig_all_combinations']], on=merge_cols)
    out['Robust'] = (out['FDR_significant'] & out['LOO_sign_always_preserved']
                     & out['LOO_all_nominal_p_lt_0.05']
                     & out['related_tool_families_FDR_sig_all_combinations']
                     & out['related_tool_families_extended_FDR_sig_all_combinations'])
    return out


robust_results = {}
for t in ROBUSTNESS_THRESHOLDS:
    elig_t = eligible_criteria(data, CRITERIA_ORDER, t)
    base_t = primary if t == PRIMARY_MIN_POSITIVE else correlate(data, elig_t)
    tag = f'min{t}of{n_obs}'
    loo_t, dep_scen_t, dep_pairs_t = robustness(elig_t, base_t, tag)
    robust_results[t] = {'eligible': elig_t, 'loo': loo_t, 'dep_scen': dep_scen_t,
                         'flagged': add_robust_flags(base_t, loo_t, dep_pairs_t)}
    la = loo_t[loo_t.FDR_family == FAMILY_A]
    print(f"   [{tag}{' primary' if t == PRIMARY_MIN_POSITIVE else ''}] family A max |Δr| "
          f"median {la['LOO_max_abs_change_r'].median():.3f}, max {la['LOO_max_abs_change_r'].max():.3f}; "
          f"largest change most often from {la['Instance_with_max_change'].value_counts().idxmax()}")
    print(dep_scen_t[['Scenario', 'N_instances', 'N_family_A_FDR_significant',
                      'N_family_B_FDR_significant']].to_string(index=False))
loo = robust_results[PRIMARY_MIN_POSITIVE]['loo']
dep_scen = robust_results[PRIMARY_MIN_POSITIVE]['dep_scen']
# ---------------------------------------------------------------------------
# 8. Primary long table with robustness flags
# ---------------------------------------------------------------------------
pos = positive_counts(data, CRITERIA_ORDER)
primary = robust_results[PRIMARY_MIN_POSITIVE]['flagged']
at = list(primary.columns).index('FDR_significant') + 1
primary.insert(at, 'N_positive_jaccard_var1', primary['Variable_1'].map(pos))
primary.insert(at + 1, 'N_positive_jaccard_var2', primary['Variable_2'].map(pos))
primary.to_csv(output_dir / 'pearson_pairwise_long.csv', index=False)

# Consistency: long table == matrices
for _, row in primary.iterrows():
    assert r_matrix.loc[row.Variable_1, row.Variable_2] == row.Pearson_r
    assert q_matrix.loc[row.Variable_2, row.Variable_1] == row.FDR_q_value

# ---------------------------------------------------------------------------
# 9. Concise analysis summary (all numbers taken from the saved tables)
# ---------------------------------------------------------------------------
def fmt(x, d=3):
    return 'NA' if pd.isna(x) else (f'{x:.2e}' if 0 < abs(x) < 1e-3 else f'{x:.{d}f}')


def md_table(frame):
    cols = list(frame.columns)
    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    for _, r in frame.iterrows():
        lines.append('| ' + ' | '.join(fmt(v) if isinstance(v, float) else str(v) for v in r) + ' |')
    return '\n'.join(lines)


lines = ['# Figure 5 correlation analysis summary', '',
         'Generated by `calculate_pairwise_pearson_correlations.py`; all values are read from the saved CSV files.', '',
         '## Design',
         f'- Observations: {n_obs} ClinGen tool instances (master table {df.shape[0]} x {df.shape[1]}, unchanged).',
         f'- Primary eligibility: strictly positive criterion-level Jaccard in >= {PRIMARY_MIN_POSITIVE}/{n_obs} '
         'instances and nonzero variance; all 14 observations (including zeros) retained. '
         'This is an exploratory threshold (not an established or prospectively prespecified cutoff); '
         '3/14 and 7/14 results are evaluated in sensitivity analyses.',
         '- Zero criterion-level Jaccard = no overlap with the reference for that criterion; it does not imply '
         'that the tool never assigned it (see `criterion_assignment_vs_jaccard_*.csv`).',
         f'- Included criteria ({len(eligible)}): {", ".join(eligible)}',
         f'- Excluded criteria ({len(excluded)}): ' + '; '.join(
             f"{r.Criterion} ({r.Exclusion_reason})" for r in qc[~qc.Eligible_primary].itertuples()),
         f'- FDR family A ({FAMILY_A}): {n_a} criterion-criterion tests, BH within family.',
         f'- FDR family B ({FAMILY_B}): {n_b} tests = {len(eligible) * 3} criterion-aggregate + 3 aggregate-aggregate, BH within family.',
         '- MJ, TJ and 3CA are used exactly as supplied (not recomputed after filtering).', '',
         '## Significant correlations (BH q < 0.05)']
leg_by_type = nonzero_variance_sensitivity.groupby('Pair_type')['FDR_significant'].agg(['size', 'sum'])
lines.append(f"- Nonzero-variance sensitivity analysis (single FDR family, 378 tests): {int(nonzero_variance_sensitivity.FDR_significant.sum())} significant "
             f"(criterion-criterion {int(leg_by_type.loc['criterion-criterion', 'sum'])}/{int(leg_by_type.loc['criterion-criterion', 'size'])}, "
             f"criterion-aggregate {int(leg_by_type.loc['criterion-aggregate', 'sum'])}/{int(leg_by_type.loc['criterion-aggregate', 'size'])}, "
             f"aggregate-aggregate {int(leg_by_type.loc['aggregate-aggregate', 'sum'])}/{int(leg_by_type.loc['aggregate-aggregate', 'size'])}).")
pa = primary[primary.FDR_family == FAMILY_A]
pb = primary[primary.FDR_family == FAMILY_B]
lines.append(f"- Primary analysis, family A: {int(pa.FDR_significant.sum())}/{n_a} "
             f"({int((pa.FDR_significant & (pa.Pearson_r > 0)).sum())} positive, "
             f"{int((pa.FDR_significant & (pa.Pearson_r < 0)).sum())} negative); robust: {int(pa.Robust.sum())}.")
lines.append(f"- Primary analysis, family B: {int(pb.FDR_significant.sum())}/{n_b} "
             f"(criterion-aggregate {int(pb[pb.Pair_type == 'criterion-aggregate'].FDR_significant.sum())}/{len(eligible) * 3}, "
             f"aggregate-aggregate {int(pb[pb.Pair_type == 'aggregate-aggregate'].FDR_significant.sum())}/3); robust: {int(pb.Robust.sum())}.")
lines += ['', '## Aggregate-aggregate correlations (family B)']
agg_tab = pb[pb.Pair_type == 'aggregate-aggregate'][['Variable_1', 'Variable_2', 'Pearson_r', 'P_value',
                                                      'FDR_q_value', 'LOO_min_r', 'LOO_max_r', 'r_without_Genebe',
                                                      'related_tool_families_prespecified_r']]
lines.append(md_table(agg_tab))
lines += ['', '## Strongest family A associations (sorted by r)']
show = ['Variable_1', 'Variable_2', 'Pearson_r', 'FDR_q_value', 'LOO_min_r', 'LOO_max_r',
        'Instance_with_max_change', 'r_without_Genebe', 'related_tool_families_prespecified_r', 'Robust']
lines.append(md_table(pa.sort_values('Pearson_r', ascending=False).head(8)[show]))
lines.append('')
lines.append(md_table(pa.sort_values('Pearson_r').head(5)[show]))
lines += ['', '## Family B criterion-aggregate associations with q < 0.05']
lines.append(md_table(pb[(pb.Pair_type == 'criterion-aggregate') & pb.FDR_significant]
                      .sort_values('Pearson_r', ascending=False)[show]))
lines += ['', '## Threshold sensitivity (family A)', md_table(thr_summary[[
    'Min_positive_jaccard', 'N_eligible_criteria', 'N_family_A_pairs_tested', 'N_family_A_FDR_significant',
    'Highest_positive_pair', 'Highest_positive_r', 'Most_negative_pair', 'Most_negative_r',
    'N_family_B_FDR_significant']]),
    '', 'Major pairs (nonzero-variance sensitivity top-10 |r| and top-5 |r| FDR-significant pairs at each robustness threshold) '
    'by threshold: `sensitivity/sensitivity_threshold_major_pairs.csv`.', '']
thr_long_df = pd.concat(thr_long, ignore_index=True)
cons = []
for t in SENSITIVITY_THRESHOLDS:
    if t == PRIMARY_MIN_POSITIVE:
        continue
    other = thr_long_df[thr_long_df.Analysis == f'threshold_min{t}of{n_obs}']
    m = primary.merge(other[['Variable_1', 'Variable_2', 'Pearson_r', 'FDR_q_value', 'FDR_significant']],
                      on=['Variable_1', 'Variable_2'], suffixes=('', '_other'))
    for fam, sub in m.groupby('FDR_family'):
        cons.append({'Comparison_threshold': f'{t}/{n_obs}', 'FDR_family': fam, 'N_shared_pairs': len(sub),
                     'Max_abs_r_difference': float((sub.Pearson_r - sub.Pearson_r_other).abs().max()),
                     'N_sig_both': int((sub.FDR_significant & sub.FDR_significant_other).sum()),
                     'N_sig_primary_only': int((sub.FDR_significant & ~sub.FDR_significant_other).sum()),
                     'N_sig_comparison_only': int((~sub.FDR_significant & sub.FDR_significant_other).sum())})
lines += ['## Consistency of primary pairs across thresholds',
          'Pearson r for a given pair is identical at every threshold (all 14 observations are used); '
          'only the eligible set, and hence the BH-FDR family size and q-values, change.', '',
          md_table(pd.DataFrame(cons)), '']
f7 = robust_results.get(7, {}).get('flagged')
if f7 is not None and PRIMARY_MIN_POSITIVE != 7:
    a7, b7 = f7[f7.FDR_family == FAMILY_A], f7[f7.FDR_family == FAMILY_B]
    lines += [f'## Sensitivity analysis at 7/{n_obs} comparison analysis',
              f"- Eligible ({len(robust_results[7]['eligible'])}): {', '.join(robust_results[7]['eligible'])}",
              f"- Family A: {int(a7.FDR_significant.sum())}/{len(a7)} q < 0.05; robust {int(a7.Robust.sum())}.",
              f"- Family B: {int(b7.FDR_significant.sum())}/{len(b7)} q < 0.05; robust {int(b7.Robust.sum())}.",
              '- Leave-one-out and tool-family outputs: `sensitivity/*_min7of14.csv`.', '']
lines += [
    '', '## Tool-instance dependence',
    '- Related families: CharGer_Local/CharGer_Online (same tool, two deployments); InterVar_2018/InterVar_2025 '
    '(two versions). Extended definition also groups AutoGVP with InterVar (AutoGVP runs a modified InterVar).',
    '- Pre-specified representatives: InterVar_2025 (most recent release), CharGer_Local (standalone local run); '
    'all representative combinations are also reported.', '',
    md_table(dep_scen[['Scenario', 'Prespecified', 'N_instances', 'N_family_A_FDR_significant',
                       'N_family_B_FDR_significant', 'MJ_vs_3CA_r', 'MJ_vs_3CA_q', 'TJ_vs_3CA_r', 'TJ_vs_3CA_q']]),
    '', 'Most similar instance pairs by 28-criterion Jaccard profile:', '',
    md_table(sim.head(5)),
    '', '## Robustness definition',
    '- `Robust` = BH q < 0.05 in its primary family, same sign and nominal p < 0.05 in all 14 leave-one-out '
    'subsamples, and BH q < 0.05 (recomputed) in every one-representative-per-family combination '
    '(primary and extended family definitions).',
    '- Leave-one-out and representative-subset analyses are internal sensitivity checks, not independent validation.',
    ]
(output_dir / 'correlation_analysis_summary.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')

print("\n8. Outputs written:")
for p in ['criterion_correlation_qc.csv', 'pearson_r_matrix.csv', 'pearson_pvalue_matrix.csv',
          'pearson_fdr_matrix.csv', 'pearson_fdr_family_matrix.csv', 'pearson_pairwise_long.csv',
          'correlation_analysis_summary.md']:
    print(f"   ✓ {output_dir / p}")
for p in sorted(sens_dir.glob('*.csv')):
    print(f"   ✓ {p}")

print("\n" + "=" * 80)
print("KEY AGGREGATE CORRELATIONS (family B)")
print("=" * 80)
print(agg_tab[['Variable_1', 'Variable_2', 'Pearson_r', 'P_value', 'FDR_q_value']].to_string(index=False))
print("\nPearson p-values are descriptive: 14 partly related tool instances, zero-inflated variables.")
print("=" * 80)
