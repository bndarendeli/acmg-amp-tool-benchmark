#!/usr/bin/env bash
# Display the shared analysis log in a visible WSL terminal.
set -eu
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
mkdir -p reproduced/session
touch reproduced/session/analysis.log
printf '\033]0;ACMG benchmark - WSL analysis monitor\007'
printf 'ACMG benchmark analysis monitor\nRepository: %s\n' "$repo_root"
printf 'Runtime selected: Ubuntu WSL2, Conda base, Python 3.12.11\n'
printf 'Commands run through tools/run_logged.sh appear below.\n'
printf 'This window displays the log; it is not an interactive analysis shell.\n\n'
tail -n 80 -F reproduced/session/analysis.log
