"""Package audited input snapshots and matching raw-output candidates."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import pandas as pd
from audit_candidate_inputs import ROOT

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-root',type=Path,required=True)
p.add_argument('--audit-root',type=Path,required=True)
a=p.parse_args()
dest=ROOT/'data'
if dest.exists():p.error('data/ already exists; refusing to overwrite the release.')
audit=a.audit_root
raw=pd.read_csv(audit/'raw-audit-01/raw_comparison.csv')
raw_inventory=json.loads((audit/'raw-audit-01/raw_inventory.json').read_text())
hashes={r['source']:r['sha256'] for r in raw_inventory}
merged_inventory=json.loads((audit/'input-audit-01/inventory.json').read_text())
hashes.update({r['source']:r['sha256'] for r in merged_inventory})
selected=[]
for (tool,dataset), group in raw[raw.canonical_full_match].groupby(['tool','dataset']):
    # Equivalent alternatives are retained in the audit, not duplicated in data/.
    ordered=sorted(group.source,key=lambda s:('notvus' in s.lower(),'/results_fix/' not in s,s))
    chosen=ordered[0]
    selected.append(dict(source=chosen,tool=tool,dataset=dataset,status='canonical_prediction_and_criteria_match',
                         equivalent_candidates=ordered[1:]))
vip_source='datasets/foxl2/results_fix/foxl2_hg19_viphl.tsv'
vip_diffs=pd.read_csv(audit/'selected-details-01/viphl_foxl2_differences.csv')
assert len(vip_diffs)==8 and vip_diffs.raw_parsed_class.eq('Unknown').all() and vip_diffs.merged_class.isna().all()
selected.append(dict(source=vip_source,tool='VIP-HL',dataset='foxl2',
                     status='criteria_match_eight_unknown_vs_missing_classifications',equivalent_candidates=[]))
manifest=[]
def copy(rel,target,kind,**extra):
    source=a.source_root/rel
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    if rel in hashes and hashes[rel]!=digest:raise ValueError(f'Source changed since audit: {rel}')
    path=ROOT/target;path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():raise ValueError(f'Destination exists: {target}')
    shutil.copy2(source,path)
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
    manifest.append(dict(source=rel,path=target,bytes=path.stat().st_size,sha256=digest,kind=kind,**extra))
for name in ('merged_results.tsv','foxl2_merged_results.tsv'):
    copy('analysis_last/data/'+name,'data/processed/'+name,'merged_predictions_with_reference')
for row in selected:
    target='data/tool_results/'+row['tool']+'/'+row['dataset']+'/'+Path(row['source']).name
    copy(row['source'],target,'matching_raw_output_candidate',tool=row['tool'],dataset=row['dataset'],
         verification=row['status'],equivalent_candidates=row['equivalent_candidates'])
metadata=ROOT/'metadata';metadata.mkdir(exist_ok=True)
(metadata/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
audit_dest=metadata/'audit';audit_dest.mkdir(exist_ok=True)
keep={'input-audit-01':['inventory.json','metric_comparison.json','metric_differences.csv','call_rate_comparison.csv','evidence_comparison.csv','criterion_jaccard_comparison.csv'],
      'reference-policy-01':['summary.json','policy_comparison.csv'],
      'selected-details-01':['summary.json','viphl_foxl2_differences.csv','merged_vs_fixed_differences.json'],
      'raw-audit-01':['raw_inventory.json','raw_comparison.csv'],
      'full-reproduction-01':['execution.json','csv_comparison.json']}
for folder,names in keep.items():
    target=audit_dest/folder;target.mkdir(exist_ok=True)
    for name in names:shutil.copy2(audit/folder/name,target/name)
pd.DataFrame(selected).to_csv(metadata/'selected_tool_results.csv',index=False)
print(f'Packaged {len(selected)} raw output candidates and two merged tables.')
print(f'Total input size: {sum(r["bytes"] for r in manifest)/1024/1024:.2f} MiB; largest: {max(r["bytes"] for r in manifest)/1024/1024:.2f} MiB.')
print('No analysis values were modified. The original source directories were not changed.')
