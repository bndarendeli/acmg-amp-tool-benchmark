"""Losslessly compress packaged raw tool outputs, retaining original hashes."""
from pathlib import Path
import gzip,hashlib,json,shutil

root=Path(__file__).resolve().parents[1]
manifest_path=root/'metadata/input_manifest.json'
manifest=json.loads(manifest_path.read_text())
for row in manifest:
    if row['kind']!='matching_raw_output_candidate' or row.get('compression')=='gzip':continue
    original=(root/row['path']).resolve()
    assert original.is_relative_to(root/'data/tool_results')
    assert hashlib.sha256(original.read_bytes()).hexdigest()==row['sha256']
    target=original.with_name(original.name+'.gz')
    if target.exists():raise ValueError(f'Destination exists: {target}')
    with original.open('rb') as src,target.open('wb') as dst:
        with gzip.GzipFile(filename='',mode='wb',fileobj=dst,mtime=0) as compressed:
            shutil.copyfileobj(src,compressed)
    with gzip.open(target,'rb') as stream:
        assert hashlib.sha256(stream.read()).hexdigest()==row['sha256']
    row.update(source_sha256=row['sha256'],source_bytes=row['bytes'],compression='gzip',
               path=target.relative_to(root).as_posix(),sha256=hashlib.sha256(target.read_bytes()).hexdigest(),bytes=target.stat().st_size)
    original.unlink()
manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Inputs now occupy {sum(x["bytes"] for x in manifest)/1024/1024:.2f} MiB; largest file {max(x["bytes"] for x in manifest)/1024/1024:.2f} MiB.')
