"""Validate archive integrity without executing archived analysis scripts."""
from pathlib import Path
import ast
import csv
import hashlib
import json
import sys

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / 'docs/source_manifest.json').read_text(encoding='utf-8'))
input_manifest_path = root / 'metadata/input_manifest.json'
input_manifest = json.loads(input_manifest_path.read_text(encoding='utf-8')) if input_manifest_path.exists() else []
errors = []
for excluded in ('analysis/fig6', 'analysis/figures', 'docs/GITHUB_PUSH_TR.md'):
    if (root / excluded).exists():
        errors.append(f'Excluded path is present: {excluded}')
for item in manifest + input_manifest:
    p = root / item['path']
    if not p.is_file():
        errors.append(f'Missing: {item["path"]}')
        continue
    if hashlib.sha256(p.read_bytes()).hexdigest() != item['sha256']:
        errors.append(f'Checksum mismatch: {item["path"]}')
    if p.stat().st_size >= 100 * 1024 * 1024:
        errors.append(f'File reaches 100 MiB: {item["path"]}')
for p in root.rglob('*.py'):
    if any(x in p.parts for x in ('.venv', 'reproduced')):
        continue
    try:
        ast.parse(p.read_text(encoding='utf-8-sig'), filename=str(p))
    except (SyntaxError, UnicodeError) as exc:
        errors.append(str(exc))
schemas = json.loads((root / 'docs/table_schema.json').read_text(encoding='utf-8'))
archived_csv_count = len(schemas)
processed_schema_path = root / 'metadata/processed_table_schema.json'
if processed_schema_path.exists():
    schemas += json.loads(processed_schema_path.read_text(encoding='utf-8'))
for schema in schemas:
    p = root / schema['path']
    if not p.is_file():
        continue
    with p.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.reader(stream, delimiter='\t' if p.suffix == '.tsv' else ',')
        columns = next(reader, [])
        count = 0
        for row in reader:
            count += 1
            if len(row) != len(columns):
                errors.append(f'CSV width mismatch: {schema["path"]}, row {count+1}')
        if columns != schema['columns'] or count != schema['rows']:
            errors.append(f'CSV schema/count mismatch: {schema["path"]}')
expected = {item['path'] for item in manifest}
actual = {p.relative_to(root).as_posix() for p in (root / 'analysis').rglob('*') if p.is_file()}
for name in sorted(actual - expected):
    errors.append(f'Untracked archive file: {name}')
if input_manifest:
    expected_inputs = {item['path'] for item in input_manifest}
    actual_inputs = {p.relative_to(root).as_posix() for p in (root / 'data').rglob('*') if p.is_file() and p.suffix.lower() != '.md'}
    for name in sorted(actual_inputs - expected_inputs):
        errors.append(f'Untracked input file: {name}')
if errors:
    print('\n'.join(errors))
    sys.exit(1)
print(f'PASS: {len(manifest)} archived files; {len(input_manifest)} input files; {archived_csv_count} archived CSVs and {len(schemas)-archived_csv_count} processed TSVs; checksums, syntax and structure verified.')
