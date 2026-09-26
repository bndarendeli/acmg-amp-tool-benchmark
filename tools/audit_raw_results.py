"""Trace candidate raw output files to normalized manuscript predictions."""
from pathlib import Path
import argparse
import contextlib
import hashlib
import io
import json
import sys
import pandas as pd
from audit_candidate_inputs import ROOT, functions_from, parse_criteria_string

sys.path.insert(0,str(ROOT/'analysis/fig2/parsers'))
from parser_factory import ParserFactory

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-root',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
args=p.parse_args()
args.output.mkdir(parents=True,exist_ok=False)
source=args.source_root
ns=functions_from(ROOT/'analysis/fig5/build_fig4_master_table.py',('CLASS_NORMALIZATION',))
def classification(x):
    if pd.isna(x) or str(x).strip()=='': return ''
    return ns['normalize_classification'](x) or str(x).strip()
def criteria(x):
    return ','.join(sorted(parse_criteria_string(x)))
mapping={'Intervar_20180118':'InterVar_2018','Intervar_20250721':'InterVar_2025','BIAS':'BIAS',
         'CharGer_local':'CharGer_Local','CharGer_online':'CharGer_Online','DiabloACMG':'DiabloACMG',
         'Exomiser':'Exomiser','Genebe':'Genebe','TAPES':'TAPES','Franklin':'Franklin','AutoGVP':'AutoGVP',
         'VIP-HL':'VIP-HL','Cancer SIGVAR':'CancerSIGVAR','CPSR':'CPSR'}
merged=pd.read_csv(source/'analysis_last/data/merged_results.tsv',sep='\t',low_memory=False)
fox=pd.read_csv(source/'analysis_last/data/foxl2_merged_results.tsv',sep='\t',low_memory=False)
cohorts={name: merged.loc[merged.Dataset==name] for name in ('clingen_28012026','HGMD_Clinvar_HL','HGMD_Clinvar_Cancer')}
cohorts['foxl2']=fox
files=[]
for folder,tool in mapping.items():
    files.extend((f,tool,False) for f in sorted((source/'Results'/folder).glob('*')) if f.suffix.lower() in ('.csv','.tsv','.txt','.vcf'))
for folder in ('results','results_fix'):
    for f in sorted((source/'datasets/foxl2'/folder).glob('*')):
        if f.suffix.lower() not in ('.csv','.tsv','.txt','.vcf'):continue
        low=f.name.lower()
        tool=next((v for token,v in [('20180118','InterVar_2018'),('20250721','InterVar_2025'),('charger_local','CharGer_Local'),
             ('charger_online','CharGer_Online'),('diablo','DiabloACMG'),('autogvp','AutoGVP'),('exomiser','Exomiser'),
             ('genebe','Genebe'),('franklin','Franklin'),('viphl','VIP-HL'),('tapes','TAPES'),('bias','BIAS')] if token in low),None)
        if tool:files.append((f,tool,True))
inventory=[]; comparisons=[]
for i,(f,tool,is_fox) in enumerate(files):
    rel=f.relative_to(source).as_posix()
    print(f'[{i+1}/{len(files)}] {rel}',flush=True)
    record=dict(source=rel,tool=tool,bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
    messages=io.StringIO()
    try:
        with contextlib.redirect_stdout(messages):
            parsed=ParserFactory.get_parser(tool,str(f)).parse(str(f))
        if parsed is None or parsed.empty:
            record.update(status='no parsed rows',parser_messages=messages.getvalue()[:2000]);inventory.append(record);continue
        record['parsed_rows']=len(parsed)
        record['duplicate_keys']=int(parsed.Variant_Key.duplicated().sum())
        parsed=parsed.drop_duplicates('Variant_Key',keep='first').set_index('Variant_Key')
        for name,sub in cohorts.items():
            if (name=='foxl2') != is_fox:continue
            key='Variant_Key_hg19' if tool in ('Franklin','VIP-HL','CancerSIGVAR') else 'Variant_Key'
            col=tool+'_Classification'; critcol=tool+'_ACMG_Criteria'
            if col not in sub:continue
            aligned=parsed.reindex(sub[key].astype(str))
            actual_class=aligned.Classification.map(classification).to_numpy()
            expected_class=sub[col].map(classification).to_numpy()
            actual_criteria=aligned.ACMG_Criteria.map(criteria).to_numpy()
            expected_criteria=sub[critcol].map(criteria).to_numpy()
            overlap=int(sub[key].isin(parsed.index).sum())
            cm=int((actual_class!=expected_class).sum()); em=int((actual_criteria!=expected_criteria).sum())
            comparisons.append(dict(source=rel,tool=tool,dataset=name,cohort_rows=len(sub),overlapping_keys=overlap,
                                    classification_mismatches=cm,criteria_mismatches=em,
                                    canonical_full_match=overlap>0 and cm==0 and em==0))
        record['status']='parsed'
    except Exception as exc:
        record.update(status='error',error=str(exc),parser_messages=messages.getvalue()[:2000])
    inventory.append(record)
    # Persist incremental results so long audits remain inspectable.
    (args.output/'raw_inventory.json').write_text(json.dumps(inventory,indent=2))
    pd.DataFrame(comparisons).to_csv(args.output/'raw_comparison.csv',index=False)
(args.output/'raw_inventory.json').write_text(json.dumps(inventory,indent=2))
result=pd.DataFrame(comparisons)
result.to_csv(args.output/'raw_comparison.csv',index=False)
matches=result[result.canonical_full_match]
print('\nFull canonical prediction/evidence matches:',flush=True)
print(matches[['source','tool','dataset','overlapping_keys']].to_string(index=False),flush=True)
print(f'Checked {len(inventory)} files; {len(matches)} file/cohort matches.',flush=True)
