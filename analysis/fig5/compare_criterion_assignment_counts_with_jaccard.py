#!/usr/bin/env python3
"""
Descriptive comparison of criterion assignment counts and criterion-level Jaccard
sparsity (Figure 5 QC companion; not used as input to any correlation).

A zero criterion-level Jaccard means that the tool's and the reference's variant
sets for that criterion do not overlap. It can arise when the tool never assigned
the criterion, or when it assigned it only to variants without the reference
criterion. This script separates those cases.

Assignment counts are derived read-only from:
  - Figure 3 rule (panel a): fig6/outputs/fig6_master_variant_tool.csv, counting
    each tool's parsed `tool_criteria` per variant row (all 28 criteria here; the
    published Figure 3 panel shows a subset).
  - Jaccard evaluation set: data/processed/merged_results_noHGMD.tsv, ClinGen rows
    with available reference criteria, exactly as in fig3/calculate_clingen_jaccard.py.
    The reconstructed intersections and unions are checked against the Figure 3
    Jaccard table, which is not modified.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

script_dir = Path(__file__).resolve().parent
analysis_dir = script_dir.parent
repo_dir = analysis_dir.parent
sys.path.insert(0, str(analysis_dir / 'fig2'))
from utils_criteria_parser import parse_criteria_string  # noqa: E402

JACCARD_PATH = analysis_dir / 'fig3' / 'jaccard_per_criterion_clingen_28012026.csv'
FIG3_COUNT_SOURCE = analysis_dir / 'fig6' / 'outputs' / 'fig6_master_variant_tool.csv'
MERGED_PATH = repo_dir / 'data' / 'processed' / 'merged_results_noHGMD.tsv'
QC_PATH = script_dir / 'criterion_correlation_qc.csv'
DATASET = 'clingen_28012026'

CRITERIA_ORDER = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]
MISSING = {'', '.', 'Unknown', 'Not_Provided', 'None', 'nan', 'NA'}


def is_available(value):
    return not pd.isna(value) and str(value).strip() not in MISSING


def criteria_sets(series):
    return [parse_criteria_string(v) if is_available(v) else set() for v in series]


print("=" * 80)
print("CRITERION ASSIGNMENT COUNTS VS CRITERION-LEVEL JACCARD SPARSITY")
print("=" * 80)

jac = pd.read_csv(JACCARD_PATH)
tools = sorted(jac['Tool'].unique())

# Figure 3 counting rule (fig6 master variant/tool table)
master = pd.read_csv(FIG3_COUNT_SOURCE, low_memory=False)
master = master[master['dataset'] == DATASET]
fig3_counts = {}
for tool, sub in master.groupby('tool'):
    sets = criteria_sets(sub['tool_criteria'])
    for c in CRITERIA_ORDER:
        fig3_counts[(tool, c)] = sum(c in s for s in sets)

# Jaccard evaluation set (same variant filter as the Figure 3 Jaccard script)
merged = pd.read_csv(MERGED_PATH, sep='\t', low_memory=False)
merged = merged[merged['Dataset'].eq(DATASET)]
merged = merged[merged['Ground_Truth_ACMG'].map(is_available)]
ref_sets = [parse_criteria_string(v) for v in merged['Ground_Truth_ACMG']]
ref_has = {c: np.array([c in s for s in ref_sets]) for c in CRITERIA_ORDER}

rows = []
for tool in tools:
    tool_sets = [parse_criteria_string(v) for v in merged[f'{tool}_ACMG_Criteria']]
    for c in CRITERIA_ORDER:
        t_has = np.array([c in s for s in tool_sets])
        inter, union = int((t_has & ref_has[c]).sum()), int((t_has | ref_has[c]).sum())
        rows.append({'Tool': tool, 'Criterion': c,
                     'Reference_assignments_in_jaccard_set': int(ref_has[c].sum()),
                     'Tool_assignments_in_jaccard_set': int(t_has.sum()),
                     'Intersection_recomputed': inter, 'Union_recomputed': union,
                     'Tool_assignment_count_fig3_rule': fig3_counts.get((tool, c), 0)})
long = pd.DataFrame(rows).merge(
    jac[['Tool', 'Criterion', 'Jaccard_Similarity', 'Intersection', 'Union']], on=['Tool', 'Criterion'])
assert len(long) == 14 * 28
assert (long['Intersection_recomputed'] == long['Intersection']).all(), "Intersection mismatch with Figure 3"
assert (long['Union_recomputed'] == long['Union']).all(), "Union mismatch with Figure 3"
print("   ✓ Recomputed intersections/unions identical to the Figure 3 Jaccard table")
long = long.drop(columns=['Intersection_recomputed', 'Union_recomputed'])


def category(r):
    if r.Jaccard_Similarity > 0:
        return 'positive_overlap'
    if r.Tool_assignments_in_jaccard_set > 0:
        return 'zero_overlap_tool_assigned_criterion'
    if r.Reference_assignments_in_jaccard_set > 0:
        return 'zero_overlap_tool_never_assigned_criterion'
    return 'no_reference_and_no_tool_assignments'


long['Jaccard_category'] = long.apply(category, axis=1)
long.to_csv(script_dir / 'criterion_assignment_vs_jaccard_sparsity.csv', index=False)

summary = []
for c in CRITERIA_ORDER:
    s = long[long['Criterion'] == c]
    cat = s['Jaccard_category'].value_counts()
    summary.append({
        'Criterion': c,
        'Reference_assignments_in_jaccard_set': int(s['Reference_assignments_in_jaccard_set'].iloc[0]),
        'N_instances_positive_jaccard': int(cat.get('positive_overlap', 0)),
        'N_instances_assigning_in_jaccard_set': int((s['Tool_assignments_in_jaccard_set'] > 0).sum()),
        'N_instances_assigning_fig3_rule': int((s['Tool_assignment_count_fig3_rule'] > 0).sum()),
        'N_zero_jaccard_but_tool_assigned': int(cat.get('zero_overlap_tool_assigned_criterion', 0)),
        'N_zero_jaccard_tool_never_assigned': int(cat.get('zero_overlap_tool_never_assigned_criterion', 0)),
        'Total_tool_assignments_in_jaccard_set': int(s['Tool_assignments_in_jaccard_set'].sum()),
        'Total_tool_assignments_fig3_rule': int(s['Tool_assignment_count_fig3_rule'].sum()),
        'Instances_zero_jaccard_but_assigned': ';'.join(
            s.loc[s['Jaccard_category'] == 'zero_overlap_tool_assigned_criterion', 'Tool']),
    })
summary = pd.DataFrame(summary)
if QC_PATH.exists():
    summary = summary.merge(pd.read_csv(QC_PATH)[['Criterion', 'Eligible_primary']], on='Criterion')
summary.to_csv(script_dir / 'criterion_assignment_vs_jaccard_summary.csv', index=False)

print(summary[['Criterion', 'Reference_assignments_in_jaccard_set', 'N_instances_positive_jaccard',
               'N_instances_assigning_in_jaccard_set', 'N_zero_jaccard_but_tool_assigned']].to_string(index=False))
print(f"\n   ✓ {script_dir / 'criterion_assignment_vs_jaccard_sparsity.csv'}")
print(f"   ✓ {script_dir / 'criterion_assignment_vs_jaccard_summary.csv'}")
print("Descriptive only: assignment counts are not substituted for Jaccard in any analysis.")
