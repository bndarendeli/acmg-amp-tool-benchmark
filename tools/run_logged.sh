#!/usr/bin/env bash
# Run from an activated analysis environment; mirror commands/output to the monitor.
set -uo pipefail
export PYTHONDONTWRITEBYTECODE=1
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root" || exit 1
mkdir -p reproduced/session
log_file="$repo_root/reproduced/session/analysis.log"
if [ "$#" -eq 0 ]; then
    echo 'Usage: bash tools/run_logged.sh COMMAND [ARGUMENTS...]' >&2
    exit 2
fi
{
    printf '\n[%s] Environment: %s\n' "$(date -Is)" "${CONDA_DEFAULT_ENV:-unactivated}"
    printf 'Command:'
    printf ' %q' "$@"
    printf '\n'
    "$@"
    command_status=$?
    printf '\nExit status: %s\n' "$command_status"
    exit "$command_status"
} 2>&1 | tee -a "$log_file"
exit "${PIPESTATUS[0]}"
