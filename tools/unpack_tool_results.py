"""Restore bundled compressed tool outputs to an ignored working directory."""
import gzip
import hashlib
import json
import os
from pathlib import Path


def unpack_tool_results(repo):
    repo = Path(repo).resolve()
    destination = repo / 'reproduced/raw_inputs'
    if os.name == 'nt' and not str(destination).startswith('\\\\?\\'):
        destination = Path('\\\\?\\' + str(destination))
    for item in json.loads((repo / 'metadata/input_manifest.json').read_text()):
        if item.get('compression') != 'gzip':
            continue
        target = (destination / item['source']).resolve()
        if not target.is_relative_to(destination.resolve()):
            raise ValueError('Input manifest path escapes the working directory')
        content = gzip.decompress((repo / item['path']).read_bytes())
        if hashlib.sha256(content).hexdigest() != item['source_sha256']:
            raise ValueError(f"Checksum mismatch: {item['path']}")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
    return destination


if __name__ == '__main__':
    print(unpack_tool_results(Path(__file__).resolve().parents[1]))
