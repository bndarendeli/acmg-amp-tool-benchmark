"""Check bundled inputs, cohort boundaries, portable paths and Python syntax."""
import ast
import csv
import hashlib
import json
from pathlib import Path
import re


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / 'metadata/input_manifest.json').read_text())
    for item in manifest:
        path = root / item['path']
        assert path.is_file(), f'Missing input: {path}'
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'], f'Changed input: {path}'
    provenance = json.loads((root / 'metadata/source_manifest.json').read_text())
    for item in provenance:
        path = root / item['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'], f'Changed supplied artifact: {path}'
    with (root / 'data/processed/merged_results_noHGMD.tsv').open(encoding='utf-8', newline='') as stream:
        counts = {}
        for row in csv.DictReader(stream, delimiter='\t'):
            counts[row['Dataset']] = counts.get(row['Dataset'], 0) + 1
    assert counts == {'clingen_28012026': 11409, 'FOXL2': 288}, counts
    scripts = list((root / 'analysis').rglob('*.py')) + list((root / 'tools').glob('*.py'))
    for path in scripts:
        source = path.read_text(encoding='utf-8-sig')
        ast.parse(source, filename=str(path))
        assert not re.search(r'/mnt/[a-z]/|[A-Z]:\\', source), f'Absolute machine path: {path}'
    for name in ['fig2', 'fig3', 'fig4', 'fig5', 'fig6']:
        assert (root / 'analysis' / name).is_dir()
    for name in ['fig2', 'fig3', 'fig4', 'fig5', 'figures']:
        assert not (root / name).exists(), f'Obsolete duplicate directory: {name}'
    assert len([m for m in manifest if m.get('compression') == 'gzip']) == 26
    print(f'PASS: {len(manifest)} inputs, {len(scripts)} Python files, cohorts {counts}.')


if __name__ == '__main__':
    main()
