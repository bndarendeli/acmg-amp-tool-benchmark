"""Compare included/excluded missing-reference policies without changing results."""
from pathlib import Path
import argparse
import json
import pandas as pd
from audit_candidate_inputs import ROOT, functions_from, parse_criteria_string

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--input',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
args=p.parse_args()
args.output.mkdir(parents=True,exist_ok=False)
df=pd.read_csv(args.input,sep='\t',low_memory=False)
df=df[df.Dataset=='clingen_28012026']
target=pd.read_csv(ROOT/'analysis/fig3/jaccard_per_criterion_clingen_28012026.csv')
available=functions_from(ROOT/'analysis/fig5/build_fig4_master_table.py')['is_criteria_available']
ref=df.Ground_Truth_ACMG.map(parse_criteria_string)
mask=df.Ground_Truth_ACMG.map(available)
cache={}; rows=[]
for t in target.itertuples(index=False):
    if t.Tool not in cache:
        cache[t.Tool]=df[t.Tool+'_ACMG_Criteria'].map(parse_criteria_string)
    a=ref.map(lambda s:t.Criterion in s); b=cache[t.Tool].map(lambda s:t.Criterion in s)
    for policy, selected in [('exclude_missing_reference',mask),('include_missing_reference',pd.Series(True,index=df.index))]:
        inter=int((a[selected]&b[selected]).sum()); union=int((a[selected]|b[selected]).sum())
        rows.append(dict(tool=t.Tool,criterion=t.Criterion,policy=policy,intersection=inter,union=union,
                         jaccard=inter/union if union else 0,archived_jaccard=t.Jaccard_Similarity,
                         match=inter==t.Intersection and union==t.Union))
result=pd.DataFrame(rows)
result.to_csv(args.output/'policy_comparison.csv',index=False)
summary={'variants':len(df),'missing_reference':int((~mask).sum()),
         'policies':result.groupby('policy').match.agg(['sum','count']).to_dict(orient='index')}
(args.output/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
