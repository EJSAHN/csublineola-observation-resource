#!/usr/bin/env bash
#SBATCH --job-name=csub-verification
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
if [[ $# -lt 3 || $# -gt 5 ]]; then
  printf '%s\n' 'Usage: sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION run_ceres.sh /absolute/repository /absolute/resource /absolute/new_output [/absolute/followup_audit] [/absolute/observation_masks]' >&2
  exit 2
fi
# Compatibility entry point; use run_slurm.sh in new submissions.
exec bash "$1/run_slurm.sh" "$@"
