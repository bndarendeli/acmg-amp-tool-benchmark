"""Inspect candidate equivalence and unresolved FOXL2 VIP-HL records."""
from pathlib import Path
import argparse,json,sys
import pandas as pd
from audit_candidate_inputs import ROOT
sys.path.insert(0,str(ROOT/'analysis/fig2/parsers'))
from viphl_parser import VIPHLParser

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-root',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
d=a.source_root/'analysis_last/data'
main=pd.read_csv(d/'merged_results.tsv',sep='\t',low_memory=False)
fixed=pd.read_csv(d/'merged_results_fixed.tsv',sep='\t',low_memory=False)
fox=pd.read_csv(d/'foxl2_merged_results.tsv',sep='\t',low_memory=False)
notes={'merged_rows':len(main),'dataset_counts':main.Dataset.fillna('<missing>').value_counts().to_dict(),
       'merged_vs_fixed_columns_equal':list(main.columns)==list(fixed.columns)}
if main.shape==fixed.shape and list(main.columns)==list(fixed.columns):
    diff=main.fillna('').ne(fixed.fillna(''))
    notes['merged_vs_fixed_different_cells_by_column']={k:int(v) for k,v in diff.sum().items() if v}
    sample=[]
    for idx,col in zip(*diff.to_numpy().nonzero()):
        sample.append({'row':int(idx),'column':main.columns[col],'merged':str(main.iloc[idx,col]),'fixed':str(fixed.iloc[idx,col])})
    (a.output/'merged_vs_fixed_differences.json').write_text(json.dumps(sample,indent=2))
raw=VIPHLParser().parse(str(a.source_root/'datasets/foxl2/results_fix/foxl2_hg19_viphl.tsv'))
notes['viphl_duplicate_keys']=int(raw.Variant_Key.duplicated().sum())
raw=raw.drop_duplicates('Variant_Key',keep='first').set_index('Variant_Key')
aligned=raw.reindex(fox.Variant_Key_hg19)
result=pd.DataFrame({'variant_hg38':fox.Variant_Key.to_numpy(),'variant_hg19':fox.Variant_Key_hg19.to_numpy(),
                     'raw_parsed_class':aligned.Classification.to_numpy(),'merged_class':fox['VIP-HL_Classification'].to_numpy()})
result=result[result.raw_parsed_class.fillna('').ne(result.merged_class.fillna(''))]
result.to_csv(a.output/'viphl_foxl2_differences.csv',index=False)
notes['viphl_label_difference_counts']=result.groupby(['raw_parsed_class','merged_class'],dropna=False).size().to_string()
(a.output/'summary.json').write_text(json.dumps(notes,indent=2))
print(json.dumps(notes,indent=2));print(result.to_string(index=False))
