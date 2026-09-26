"""Run supported original plot scripts on an isolated copy of archived data."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--figure', choices=('3', '4'), required=True)
args = parser.parse_args()
scripts = {'3': 'fig3/create_figure3_criterion_counts_and_jaccard.py',
           '4': 'fig4/create_figure4_dotplot.py'}
run = root / 'reproduced' / ('figure-' + args.figure)
if run.exists():
    parser.error(f'Run directory already exists: {run}. Move it before starting a fresh run.')
working = run / 'analysis'
shutil.copytree(root / 'analysis', working)
script = working / scripts[args.figure]
env = dict(os.environ, MPLBACKEND='Agg', PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1',
           MPLCONFIGDIR=str(run / '.mplconfig'))
result = subprocess.run([sys.executable, '-X', 'utf8', str(script)], cwd=script.parent, env=env)
if result.returncode:
    raise SystemExit(result.returncode)
print(f'Regenerated output is inside {working}; archived originals were not changed.')
