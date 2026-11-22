#!/usr/bin/env bash
#SBATCH --job-name=csub-verification
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
if [[ $# -lt 3 || $# -gt 5 ]]; then
  printf '%s\n' 'Usage: sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION run_slurm.sh /absolute/repository /absolute/resource /absolute/new_output [/absolute/followup_audit] [/absolute/observation_masks]' >&2
  exit 2
fi
repo_dir="$1"
data_dir="$2"
output_dir="$3"
# Slurm copies this script to its spool, so use the supplied repository path.
args=("$repo_dir/run.py" --data "$data_dir" --out "$output_dir")
if [[ $# -ge 4 && -n "$4" ]]; then args+=(--followup "$4"); fi
if [[ $# -eq 5 && -n "$5" ]]; then args+=(--masks "$5"); fi
export PYTHONDONTWRITEBYTECODE=1
python3 "${args[@]}"
