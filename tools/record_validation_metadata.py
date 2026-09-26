"""Record the validated local analysis environment and packaged-input catalog."""
from pathlib import Path
import importlib.metadata
import json
import platform
import shutil
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[1]
packages=('pandas','numpy','matplotlib','seaborn','scikit-learn','adjustText','openpyxl')
runtime={'python':platform.python_version(),'platform':platform.platform(),
         'packages':{p:importlib.metadata.version(p) for p in packages}}
(root/'metadata/runtime_versions.json').write_text(json.dumps(runtime,indent=2)+'\n')
manifest=json.loads((root/'metadata/input_manifest.json').read_text())
policy=pd.read_csv(root/'metadata/audit/reference-policy-01/policy_comparison.csv')
excluded=policy[policy.policy=='exclude_missing_reference']
changed_values=int((~np.isclose(excluded.jaccard,excluded.archived_jaccard,rtol=0,atol=1e-12)).sum())
rawcheck=root/'reproduced/packaged-raw-rerun-01/comparison.csv'
if not rawcheck.exists():raise SystemExit('Packaged raw reparse is still incomplete.')
raw=pd.read_csv(rawcheck)
assert len(raw)==52 and not raw[['classification_mismatches','criteria_mismatches']].to_numpy().any()
shutil.copy2(rawcheck,root/'metadata/audit/packaged_raw_reparse.csv')
size=sum(x['bytes'] for x in manifest)/1024/1024
validation=f'''# Validation results

Validated on 2026-09-12 using Ubuntu WSL2, Conda base, Python {runtime['python']}.

- 125 archived analysis files and 33 archived CSV tables; original scientific values retained.
- 54 added input files: 52 gzip-compressed raw-output candidates and two merged TSVs.
- Added inputs total {size:.2f} MiB; largest stored input is {max(x['bytes'] for x in manifest)/1024/1024:.2f} MiB.
- Gzip decompression was checked against original source SHA-256 hashes during compression.
- All 52 populated tool/cohort combinations were reparsed from the packaged files and matched normalized five-tier classifications and canonical criteria. Two all-no-call combinations were separately verified. Eight FOXL2 VIP-HL Unknown/missing encoding differences remain documented.
- All 594 Figure 4 metric cells and 54 Figure 2 call counts matched the selected inputs.
- Twelve analysis/plotting scripts executed successfully in an isolated copy. Rewritten CSV outputs matched the archive; images are not claimed to be byte-identical.
- ClinGen Figure 3: excluding 257 missing-reference variants changes 109 intersection/union count pairs and {changed_values} Jaccard values. The archived table instead matches inclusion of those rows. This inconsistency is preserved for author review.
- No original source files were edited; no remote repository, commit or push was created.

Run `python tools/validate_repository.py` for current integrity checks. See `INPUT_VERIFICATION.md` for comparison scope and limitations, and `metadata/audit/` for machine-readable evidence.
'''
(root/'docs/VALIDATION.md').write_text(validation)
catalog=['# Processed and raw input catalog','','Two processed TSVs contain reference data and predictions. Counts below describe physical records, not automatically the evaluated population.','']
schema=[]
for row in manifest:
    if row['kind']!='merged_predictions_with_reference':continue
    frame=pd.read_csv(root/row['path'],sep='\t',low_memory=False)
    schema.append(dict(path=row['path'],rows=len(frame),columns=list(frame.columns)))
    catalog += [f'## [{Path(row["path"]).name}](../{row["path"]})','',f'Rows: {len(frame):,}; columns: {len(frame.columns)}.','',
                'Columns: '+', '.join('`'+c+'`' for c in frame.columns),'']
catalog+=['## Raw tool outputs','','| Tool | Cohort | File | Verification |','| --- | --- | --- | --- |']
for row in manifest:
    if row['kind']=='matching_raw_output_candidate':
        catalog.append(f'| {row["tool"]} | {row["dataset"]} | [{Path(row["path"]).name}](../{row["path"]}) | {row["verification"]} |')
(root/'docs/INPUT_CATALOG.md').write_text('\n'.join(catalog)+'\n')
(root/'metadata/processed_table_schema.json').write_text(json.dumps(schema,indent=2)+'\n')
print(f'Metadata recorded. Missing-reference policy changes {changed_values} Jaccard values; raw reparse passed for all 52 populated combinations.')
