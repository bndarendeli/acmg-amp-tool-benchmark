"""Reparse bundled raw outputs and compare normalized predictions to snapshots.

Reference fields and coordinate mappings come from the bundled merged tables.
No archived scientific values are overwritten.
"""
from pathlib import Path
import argparse,gzip,json,shutil,sys,tempfile
import numpy as np
import pandas as pd
from audit_candidate_inputs import ROOT,functions_from,parse_criteria_string
sys.path.insert(0,str(ROOT/'analysis/fig2/parsers'))
from parser_factory import ParserFactory

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
manifest=json.loads((ROOT/'metadata/input_manifest.json').read_text())
ns=functions_from(ROOT/'analysis/fig5/build_fig4_master_table.py',('CLASS_NORMALIZATION',))
def norm(x):return ns['normalize_classification'](x) or '<no-call>'
def crit(x):return ','.join(sorted(parse_criteria_string(x)))
main=pd.read_csv(ROOT/'data/processed/merged_results.tsv',sep='\t',low_memory=False)
fox=pd.read_csv(ROOT/'data/processed/foxl2_merged_results.tsv',sep='\t',low_memory=False)
cohorts={name:main[main.Dataset==name].copy() for name in ('clingen_28012026','HGMD_Clinvar_HL','HGMD_Clinvar_Cancer')}
cohorts['foxl2']=fox.copy()
rebuilt={k:v[[c for c in v.columns if not c.endswith(('_Classification','_ACMG_Criteria')) or c.startswith('Ground_Truth_')]].copy() for k,v in cohorts.items()}
checks=[]
for row in manifest:
    if row['kind']!='matching_raw_output_candidate':continue
    tool=row['tool'];name=row['dataset'];source=ROOT/row['path']
    print(f'Reparsing {tool} / {name}',flush=True)
    with tempfile.TemporaryDirectory(dir=a.output) as tmp:
        original=Path(tmp)/Path(row['source']).name
        if row.get('compression')=='gzip':
            with gzip.open(source,'rb') as src,original.open('wb') as dst:shutil.copyfileobj(src,dst)
        else:shutil.copy2(source,original)
        parsed=ParserFactory.get_parser(tool,str(original)).parse(str(original))
    if parsed is None or parsed.empty:raise ValueError(f'No parsed records: {row["path"]}')
    parsed=parsed.drop_duplicates('Variant_Key',keep='first').set_index('Variant_Key')
    ref=cohorts[name]
    key='Variant_Key_hg19' if tool in ('Franklin','VIP-HL','CancerSIGVAR') else 'Variant_Key'
    aligned=parsed.reindex(ref[key].astype(str))
    col=tool+'_Classification';criteria_col=tool+'_ACMG_Criteria'
    rebuilt[name][col]=aligned.Classification.to_numpy()
    rebuilt[name][criteria_col]=aligned.ACMG_Criteria.to_numpy()
    cm=int(np.count_nonzero(ref[col].map(norm).to_numpy()!=aligned.Classification.map(norm).to_numpy()))
    em=int(np.count_nonzero(ref[criteria_col].map(crit).to_numpy()!=aligned.ACMG_Criteria.map(crit).to_numpy()))
    checks.append(dict(tool=tool,dataset=name,rows=len(ref),classification_mismatches=cm,criteria_mismatches=em,
                       comparison='five-tier classification; no-call encodings collapsed; canonical criterion sets'))
# Two tool/cohort combinations have no predictions in the benchmark.
for name,tool in [('HGMD_Clinvar_HL','CancerSIGVAR'),('HGMD_Clinvar_Cancer','VIP-HL')]:
    ref=cohorts[name];col=tool+'_Classification';cc=tool+'_ACMG_Criteria'
    assert ref[col].map(norm).eq('<no-call>').all() and ref[cc].map(crit).eq('').all()
    rebuilt[name][col]=np.nan;rebuilt[name][cc]=np.nan
for name,frame in rebuilt.items():frame.to_csv(a.output/(name+'_reparsed.tsv'),sep='\t',index=False)
pd.DataFrame(checks).to_csv(a.output/'comparison.csv',index=False)
failures=[r for r in checks if r['classification_mismatches'] or r['criteria_mismatches']]
print(f'{len(checks)} populated tool/cohort combinations checked; {len(failures)} mismatches. Two all-no-call combinations verified separately.',flush=True)
print('This verifies normalized content, not byte identity or original run history. Eight FOXL2 VIP-HL no-call encodings differ as documented.',flush=True)
if failures:raise SystemExit(1)
