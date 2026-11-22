#!/usr/bin/env bash
#SBATCH --job-name=csub-resource-replay
#SBATCH --partition=ceres
#SBATCH --account=ahn_project_1_ml
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
if [[ $# -lt 3 || $# -gt 5 ]]; then
  printf '%s\n' 'Usage: sbatch run_ceres.sh /absolute/repository /absolute/release_candidate /absolute/new_output [/absolute/completed_SD13_run] [/absolute/recovered_mask_data]' >&2
  exit 2
fi
repo_dir="$1"
data_dir="$2"
output_dir="$3"
# Pass an absolute repository directory because Slurm copies this script to its spool.
args=("$repo_dir/run.py" --data "$data_dir" --out "$output_dir")
if [[ $# -ge 4 && -n "$4" ]]; then args+=(--sd13 "$4"); fi
if [[ $# -eq 5 ]]; then args+=(--maskdata "$5"); fi
export PYTHONDONTWRITEBYTECODE=1
python3 "${args[@]}"
