"""Run the manuscript analyses in a separate, self-contained output copy."""
import argparse
import datetime
import os
from pathlib import Path
import shutil
import subprocess
import sys

STEPS = [
    'fig6/build_fig4_master_table.py',
    'fig6/build_capability_matrix.py',
    'fig6/perform_capability_adjusted_analysis.py',
    'fig6/perform_foxl2_capability_adjusted_analysis.py',
    'fig6/build_evidence_distributions.py',
    'fig6/calculate_total_criterion_jaccard.py',
    'fig6/calculate_foxl2_total_criterion_jaccard.py',
    'fig2/recalculate_call_rates.py',
    'fig2/create_figure2_comprehensive.py',
    'fig2/create_supp_intervar_comparison.py',
    'fig3/calculate_clingen_jaccard.py',
    'fig3/calculate_foxl2_jaccard.py',
    'fig3/create_figure3_criterion_counts_and_jaccard.py',
    'fig4/scripts/reconstruct_classification_metrics.py',
    'fig4/create_figure4_reconstructed.py',
    'fig4/create_figure4_supplementary.py',
    'fig5/build_correlation_master_table.py',
    'fig5/calculate_pairwise_pearson_correlations.py',
    'fig5/compare_criterion_assignment_counts_with_jaccard.py',
    'fig5/plot_clustered_pairwise_pearson_heatmap.py',
    'fig6/generate_supplementary_figure_s4.py',
    'fig6/generate_figure6.py',
]


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='New directory for the analysis copy and logs')
    parser.add_argument('--from-raw', action='store_true', help='Rebuild merged predictions from bundled tool outputs first')
    args = parser.parse_args()
    output = (args.output or root / 'reproduced' / datetime.datetime.now().strftime('%Y%m%d-%H%M%S')).resolve()
    if os.name == 'nt' and not str(output).startswith('\\\\?\\'):
        output = Path('\\\\?\\' + str(output))
    if output.exists():
        parser.error('Output directory already exists; choose a new directory.')
    for name in ['analysis', 'data', 'tools', 'metadata']:
        shutil.copytree(root / name, output / name, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    steps = list(STEPS)
    if args.from_raw:
        steps.insert(0, 'fig2/parse_all_results.py')
    logdir = output / 'logs'
    logdir.mkdir()
    env = dict(os.environ, MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1',
               MPLCONFIGDIR=str(output / '.mplconfig'))
    print(f'Output: {output}', flush=True)
    for index, rel in enumerate(steps, 1):
        script = output / 'analysis' / rel
        logfile = logdir / f'{index:02d}-{script.stem}.log'
        print(f'[{index}/{len(steps)}] {rel}', flush=True)
        with logfile.open('w', encoding='utf-8') as log:
            result = subprocess.run([sys.executable, str(script)], cwd=script.parent,
                                    env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            print(logfile.read_text(encoding='utf-8')[-5000:])
            raise SystemExit(f'Failed; see {logfile}')
    print(f'Completed {len(steps)} selected manuscript-figure scripts.', flush=True)


if __name__ == '__main__':
    main()
