# Local WSL execution

The author workstation has an Ubuntu WSL2 distribution with these Conda environments:

- `base`: Python 3.12.11; pandas, numpy, matplotlib, seaborn, scikit-learn, adjustText and openpyxl imported successfully.
- `acmg-benchmarking`: an import check stopped at missing `adjustText`.

The current local execution choice is `base`. No packages were installed or modified. This is a local environment observation, not a portable environment lock.

After activating the selected environment, run commands from the repository root through the logging helper:

```bash
conda activate base
bash tools/run_logged.sh python tools/validate_repository.py
```

In a separate WSL terminal, display the command/output log:

```bash
bash tools/monitor_wsl.sh
```

The monitor tails `reproduced/session/analysis.log` and displays only commands run through the logging helper. The log is excluded from Git. Closing the monitor does not stop an analysis running separately.
