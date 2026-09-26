"""Reproduce Figures 2-5 on isolated copies of selected merged inputs.

Archived Jaccard input tables are used for Figure 3 plotting; this does not
resolve the documented missing-reference policy discrepancy.
"""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--merged',type=Path,default=ROOT/'data/processed/merged_results.tsv')
    p.add_argument('--foxl2',type=Path,default=ROOT/'data/processed/foxl2_merged_results.tsv')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    out=args.output.resolve()
    if out.exists():p.error('Output exists; choose a fresh directory.')
    work=out/'analysis'
    shutil.copytree(ROOT/'analysis',work)
    inputs=out/'inputs';inputs.mkdir()
    shutil.copy2(args.merged,inputs/'merged_results.tsv')
    shutil.copy2(args.foxl2,inputs/'foxl2_merged_results.tsv')
    def edit(relative,replacements):
        f=work/relative;s=f.read_text(encoding='utf-8')
        for old,new in replacements:
            if old not in s:raise ValueError(f'Missing patch anchor in {relative}: {old}')
            s=s.replace(old,new)
        f.write_text(s,encoding='utf-8')
    edit('fig4/scripts/calculate_metrics_from_variants.py',[("Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last/data')",f'Path({str(inputs)!r})')])
    for name in ('create_figure2_comprehensive.py','create_supp_intervar_comparison.py'):
        replacements=[("'/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last/data/merged_results.tsv'",repr(str(inputs/'merged_results.tsv')))]
        if name=='create_figure2_comprehensive.py':
            replacements.append(("'/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last/results/corrected_metrics/corrected_call_rates.csv'",repr(str(work/'fig2/corrected_call_rates.csv'))))
        edit('fig2/'+name,replacements)
    edit('fig5/build_fig4_master_table.py',[
        ("BASE_DIR = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep')",f'BASE_DIR = Path({str(out)!r})'),
        ("DATA_FILE = BASE_DIR / 'analysis_last' / 'data' / 'merged_results.tsv'",f'DATA_FILE = Path({str(inputs/"merged_results.tsv")!r})'),
        ("FOXL2_FILE = BASE_DIR / 'analysis_last' / 'foxl2_merged_results.tsv'",f'FOXL2_FILE = Path({str(inputs/"foxl2_merged_results.tsv")!r})')])
    for f in (work/'fig5').glob('*.py'):
        s=f.read_text(encoding='utf-8')
        if " / 'fig1'" in s:f.write_text(s.replace(" / 'fig1'"," / 'fig2'"),encoding='utf-8')
    scripts=['fig4/scripts/calculate_metrics_from_variants.py','fig5/build_fig4_master_table.py',
             'fig5/perform_capability_adjusted_analysis.py','fig5/perform_foxl2_capability_adjusted_analysis.py',
             'fig5/generate_foxl2_criterion_decomposition.py','fig5/generate_clingen_supplementary_tables.py',
             'fig2/create_figure2_comprehensive.py','fig2/create_supp_intervar_comparison.py',
             'fig3/create_figure3_criterion_counts_and_jaccard.py','fig4/create_figure4_dotplot.py',
             'fig5/generate_foxl2_supplementary_figure.py','fig5/generate_figure4_refined.py']
    env=dict(os.environ,MPLBACKEND='Agg',MPLCONFIGDIR=str(out/'.mplconfig'),PYTHONDONTWRITEBYTECODE='1')
    logs=out/'logs';logs.mkdir()
    statuses=[]
    for script in scripts:
        print(f'Running {script}',flush=True)
        f=work/script
        with (logs/(f.stem+'.log')).open('w') as log:
            proc=subprocess.run([sys.executable,str(f)],cwd=f.parent,env=env,stdout=log,stderr=subprocess.STDOUT)
        statuses.append(dict(script=script,exit_code=proc.returncode))
        print(f'  Exit status {proc.returncode}; log: {logs/(f.stem+".log")}',flush=True)
        if proc.returncode:
            print((logs/(f.stem+'.log')).read_text()[-3000:],flush=True)
            break
    # Compare all CSV outputs to their frozen counterparts, independent of byte serialization.
    comparisons=[]
    for original in sorted((ROOT/'analysis').rglob('*.csv')):
        rel=original.relative_to(ROOT/'analysis');generated=work/rel
        if not generated.exists():continue
        a=pd.read_csv(original);b=pd.read_csv(generated)
        try:
            pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=False,atol=1e-10,rtol=1e-9)
            match=True;error=''
        except AssertionError as exc:
            match=False;error=str(exc)[:1200]
        comparisons.append(dict(path=rel.as_posix(),equivalent=match,detail=error))
    (out/'execution.json').write_text(json.dumps(statuses,indent=2))
    (out/'csv_comparison.json').write_text(json.dumps(comparisons,indent=2))
    print('Changed CSV comparisons:',flush=True)
    print(json.dumps([c for c in comparisons if not c['equivalent']],indent=2),flush=True)
    print('Note: unchanged copied CSVs are not evidence of recalculation; consult execution.json and per-script logs.',flush=True)
    if any(s['exit_code'] for s in statuses):raise SystemExit(1)

if __name__=='__main__':main()
