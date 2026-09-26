"""Compare external candidate inputs with frozen manuscript analysis tables.

Reads external files only. Writes audit products to an explicit output directory.
"""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis/fig2'))
from utils_criteria_parser import parse_criteria_string


def functions_from(path, constants=()):
    """Load function definitions/literal constants without top-level I/O."""
    tree = ast.parse(path.read_text(encoding='utf-8-sig'))
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) or
             isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id in constants for t in n.targets)]
    namespace = dict(pd=pd, np=np, parse_criteria_string=parse_criteria_string)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.source_root.resolve()
    out = args.output.resolve()
    if out.exists():
        parser.error('Output already exists; use a fresh directory.')
    out.mkdir(parents=True)
    metric = functions_from(ROOT/'analysis/fig4/scripts/calculate_metrics_from_variants.py', ('THREE_CLASS_MAP','VALID_ACMG_CLASSES'))
    evidence = functions_from(ROOT/'analysis/fig5/build_fig4_master_table.py', ('CLASS_NORMALIZATION',))
    norm = metric['normalize_classification']
    norm5 = evidence['normalize_classification']
    available = evidence['is_criteria_available']
    targets = {'clingen_28012026': ('ClinGen', 11409, 'ClinGen_corrected_metrics.csv'),
               'HGMD_Clinvar_HL': ('HGMD+ClinVar HL', 1948, 'HGMD+ClinVar_HL_corrected_metrics.csv'),
               'HGMD_Clinvar_Cancer': ('HGMD+ClinVar Cancer', 1788, 'HGMD+ClinVar_Cancer_corrected_metrics.csv'),
               'foxl2': ('FOXL2', 288, 'FOXL2_corrected_metrics.csv')}
    master = pd.read_csv(ROOT/'analysis/fig5/outputs/fig4_master_variant_tool.csv', keep_default_na=False, dtype=str)
    call_target = pd.read_csv(ROOT/'analysis/fig2/corrected_call_rates.csv')
    candidates = sorted(set(source.glob('analysis_last/data/*merged*.tsv')) |
                        set(source.glob('analysis_2/results/*/*merged*.tsv')) |
                        set(source.glob('datasets/foxl2/analysis/results/*merged*.tsv')) |
                        set(source.glob('analysis_last/*merged*.tsv')))
    inventory, comparison, diffs, call_diffs, evidence_rows, jaccard_rows = [], [], [], [], [], []
    for p in candidates:
        rel = p.relative_to(source).as_posix()
        print(f'Inspecting {rel}', flush=True)
        df = pd.read_csv(p, sep='\t', low_memory=False)
        if 'Dataset' not in df and 'foxl2' in p.name.lower():
            df['Dataset'] = 'foxl2'
        entry = dict(source=rel, sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size,
                     rows=len(df), columns=list(df.columns), datasets=df['Dataset'].value_counts().to_dict() if 'Dataset' in df else {})
        inventory.append(entry)
        if 'Dataset' not in df or 'Ground_Truth_Classification' not in df:
            entry['status'] = 'missing required dataset/reference columns'
            continue
        for dataset, (label, expected_n, filename) in targets.items():
            sub = df.loc[df.Dataset == dataset].copy()
            if sub.empty:
                continue
            gt = sub.Ground_Truth_Classification.map(norm)
            target = pd.read_csv(ROOT/'analysis/fig4/data/corrected'/filename)
            tools = list(target.Tool)
            record = dict(source=rel, dataset=dataset, rows=len(sub), expected_rows=expected_n,
                          duplicate_variant_keys=int(sub.Variant_Key.duplicated().sum()) if 'Variant_Key' in sub else None,
                          invalid_reference=int(gt.isna().sum()), reference_counts=sub.Ground_Truth_Classification.value_counts().to_dict())
            if gt.isna().any():
                record['status'] = 'invalid reference labels; metric comparison not performed'
                comparison.append(record)
                continue
            missing_tools = []
            mismatches, tested, calculated = 0, 0, []
            for tool in tools:
                col = tool+'_Classification'
                if col not in sub:
                    missing_tools.append(tool)
                    continue
                result = metric['calculate_metrics_for_tool'](gt, sub[col].map(norm), len(sub))
                calculated.append(dict(Tool=tool, **result))
                saved = target.set_index('Tool').loc[tool]
                for key, val in result.items():
                    expected = saved[key]
                    # Archived metrics are serialized to four decimal places.
                    rounded = float(f'{val:.4f}')
                    same = pd.isna(rounded) and pd.isna(expected) or np.isclose(rounded, expected, atol=1e-10, rtol=0)
                    tested += 1
                    if not same:
                        mismatches += 1
                        diffs.append(dict(source=rel,dataset=dataset,tool=tool,metric=key,recalculated=val,archived=expected))
                ctarget = call_target[(call_target.Dataset == label) & (call_target.Tool == tool)]
                if not ctarget.empty:
                    # Reproduce archived Figure 2 script's actual non-null/non-Unknown rule.
                    raw = sub[col]
                    count = int((raw.notna() & raw.ne('') & raw.ne('Unknown')).sum())
                    strict_count = int(sub[col].map(norm).notna().sum())
                    archived = int(ctarget.iloc[0].Classified)
                    call_diffs.append(dict(source=rel,dataset=dataset,tool=tool,archived_count=archived,
                                           script_count=count,strict_count=strict_count,
                                           script_match=count==archived,strict_match=strict_count==archived))
            record.update(metric_cells=tested,metric_mismatches=mismatches,missing_tools=missing_tools,
                          status='match' if mismatches==0 and not missing_tools and len(sub)==expected_n else 'different')
            comparison.append(record)
            safe = rel.replace('/','__').replace('.tsv','')+'__'+dataset
            pd.DataFrame(calculated).to_csv(out/(safe+'_metrics.csv'),index=False)
            print(f'  {dataset}: n={len(sub)}, metric mismatches={mismatches}/{tested}, missing tools={missing_tools}',flush=True)
            if dataset not in ('clingen_28012026','foxl2'):
                continue
            ref5 = sub.Ground_Truth_Classification.map(norm5)
            refsets = sub.get('Ground_Truth_ACMG',pd.Series('',index=sub.index)).map(lambda x: ','.join(sorted(parse_criteria_string(x))) if available(x) else '')
            for tool in sorted(master.tool.unique()):
                col = tool+'_Classification'
                if col not in sub:
                    continue
                pred = sub[col].map(norm5)
                mask = ref5.notna() & pred.notna()
                crit = sub.get(tool+'_ACMG_Criteria',pd.Series('',index=sub.index)).map(lambda x: ','.join(sorted(parse_criteria_string(x))) if available(x) else '')
                actual = pd.DataFrame({'variant_id':sub.loc[mask,'Variant_Key'].astype(str), 'reference_class':ref5[mask],
                                       'tool_class':pred[mask], 'reference_criteria':refsets[mask], 'tool_criteria':crit[mask]})
                saved = master[(master.dataset==dataset)&(master.tool==tool)][actual.columns]
                if actual.variant_id.duplicated().any() or saved.variant_id.duplicated().any():
                    evidence_rows.append(dict(source=rel,dataset=dataset,tool=tool,status='duplicate keys'))
                    continue
                a=actual.set_index('variant_id'); b=saved.set_index('variant_id')
                common=a.index.intersection(b.index)
                cell_diff=int(a.loc[common].ne(b.loc[common]).sum().sum())
                evidence_rows.append(dict(source=rel,dataset=dataset,tool=tool,actual_rows=len(a),archived_rows=len(b),
                                          missing_keys=len(b.index.difference(a.index)),extra_keys=len(a.index.difference(b.index)),
                                          mismatched_cells=cell_diff,status='match' if a.index.symmetric_difference(b.index).empty and cell_diff==0 else 'different'))
            jfile = ROOT/'analysis/fig3'/('jaccard_per_criterion_clingen_28012026.csv' if dataset=='clingen_28012026' else 'jaccard_per_criterion_foxl2.csv')
            jtarget = pd.read_csv(jfile)
            if {'Tool','Criterion','Intersection','Union'}.issubset(jtarget.columns):
                validref = sub.get('Ground_Truth_ACMG',pd.Series('',index=sub.index)).map(available)
                rs = sub.loc[validref,'Ground_Truth_ACMG'].map(parse_criteria_string)
                cache={}
                for row in jtarget.itertuples(index=False):
                    col=row.Tool+'_ACMG_Criteria'
                    if col not in sub:
                        continue
                    if row.Tool not in cache:
                        cache[row.Tool]=sub.loc[validref,col].map(parse_criteria_string)
                    left=rs.map(lambda s: row.Criterion in s)
                    right=cache[row.Tool].map(lambda s: row.Criterion in s)
                    inter=int((left&right).sum()); union=int((left|right).sum())
                    jaccard_rows.append(dict(source=rel,dataset=dataset,tool=row.Tool,criterion=row.Criterion,
                                             intersection=inter,union=union,archived_intersection=row.Intersection,archived_union=row.Union,
                                             match=inter==row.Intersection and union==row.Union))
    (out/'inventory.json').write_text(json.dumps(inventory,indent=2,default=int)+'\n')
    (out/'metric_comparison.json').write_text(json.dumps(comparison,indent=2,default=int)+'\n')
    for name, rows in [('metric_differences',diffs),('call_rate_comparison',call_diffs),('evidence_comparison',evidence_rows),('criterion_jaccard_comparison',jaccard_rows)]:
        pd.DataFrame(rows).to_csv(out/(name+'.csv'),index=False)
    print(f'Audit complete: {len(candidates)} candidate files. Output: {out}',flush=True)


if __name__ == '__main__':
    main()
