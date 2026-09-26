"""Audit the historical ground-truth join used for ClinGen Jaccard counts."""
from pathlib import Path
import argparse
import hashlib
import json
import pandas as pd
from audit_candidate_inputs import ROOT, functions_from, parse_criteria_string

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=False)
ns=functions_from(ROOT/'analysis/fig3/generate_clingen_jaccard_heatmaps.py')
ns['Path']=Path
target=pd.read_csv(ROOT/'analysis/fig3/jaccard_per_criterion_clingen_28012026.csv')
summary=[]
for rel in ('analysis_2/results/all_variants/merged_results.tsv','analysis_last/data/merged_results.tsv'):
    df=pd.read_csv(args.source_root/rel,sep='\t',low_memory=False)
    for gtpath in ('analysis_2/ground_truth','analysis_last/ground_truth'):
        joined=ns['load_dataset_with_tools']('clingen_28012026',df,args.source_root/gtpath)
        valid=joined.Ground_Truth_ACMG.map(ns['is_criteria_available'])
        ref=joined.loc[valid,'Ground_Truth_ACMG'].map(parse_criteria_string)
        cache={}; details=[]
        for row in target.itertuples(index=False):
            if row.Tool not in cache:
                cache[row.Tool]=joined.loc[valid,row.Tool+'_ACMG_Criteria'].map(parse_criteria_string)
            a=ref.map(lambda s:row.Criterion in s); b=cache[row.Tool].map(lambda s:row.Criterion in s)
            inter=int((a&b).sum()); union=int((a|b).sum())
            details.append(dict(tool=row.Tool,criterion=row.Criterion,intersection=inter,union=union,
                                archived_intersection=row.Intersection,archived_union=row.Union,
                                match=inter==row.Intersection and union==row.Union))
        result=pd.DataFrame(details)
        record=dict(merged_source=rel,ground_truth_source=gtpath+'/clingen_28012026_ground_truth.xlsx',
                    joined_rows=len(joined),unique_variants=joined.Variant_Key.nunique(),
                    available_reference_rows=int(valid.sum()),matching_pairs=int(result.match.sum()),total_pairs=len(result))
        summary.append(record)
        result.to_csv(args.output/(rel.replace('/','__')+'__'+gtpath.replace('/','__')+'.csv'),index=False)
        print(record,flush=True)
(args.output/'summary.json').write_text(json.dumps(summary,indent=2))
