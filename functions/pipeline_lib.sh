#!/usr/bin/env bash
# Shared paths, logging and HPC submission helpers for the numbered
# drivers 01_-12_ in the repository root.

set -euo pipefail

BRAINOMICS_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BRAINOMICS_REPO_ROOT="$(cd "$BRAINOMICS_LIB_DIR/.." && pwd)"

# Optional site-specific configuration for the Slurm helpers (hpc/local.env is
# git-ignored; hpc/local.env.example documents the variables).
if [ -f "$BRAINOMICS_REPO_ROOT/hpc/local.env" ]; then
  # shellcheck disable=SC1091
  source "$BRAINOMICS_REPO_ROOT/hpc/local.env"
fi
# The same R/Python prefix is used by local drivers and Slurm jobs.
source "$BRAINOMICS_REPO_ROOT/environment/activate.sh"
BRAINOMICS_DATA_ROOT="${BRAINOMICS_DATA_ROOT:-$(cd "$BRAINOMICS_REPO_ROOT/../.." && pwd)/data/BrainOmicsData}"
BRAINOMICS_RUN_ROOT="${BRAINOMICS_RUN_ROOT:-$BRAINOMICS_REPO_ROOT/results/run_root}"
BRAINOMICS_RESULTS_DIR="${BRAINOMICS_RESULTS_DIR:-$BRAINOMICS_RUN_ROOT/integration}"
BRAINOMICS_EXECUTOR="${BRAINOMICS_EXECUTOR:-local}"
BRAINOMICS_RSCRIPT="${BRAINOMICS_RSCRIPT:-Rscript}"
BRAINOMICS_PYTHON="${BRAINOMICS_PYTHON:-python3}"
BRAINOMICS_SCVI_PYTHON="${BRAINOMICS_SCVI_PYTHON:-python3}"
export BRAINOMICS_REPO_ROOT BRAINOMICS_DATA_ROOT BRAINOMICS_RUN_ROOT
export BRAINOMICS_RESULTS_DIR BRAINOMICS_EXECUTOR BRAINOMICS_RSCRIPT BRAINOMICS_SCVI_PYTHON
export BRAINOMICS_PYTHON
export BRAINOMICS_HPC_ROOT BRAINOMICS_HPC_PARTITION BRAINOMICS_HPC_GPU_PARTITION
export BRAINOMICS_HPC_LOG_DIR BRAINOMICS_HPC_HOST BRAINOMICS_HPC_PORT BRAINOMICS_HPC_DATA_ROOT

brainomics_log() {
  printf '[%s] [%s] %s\n' \
    "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    "${BRAINOMICS_STAGE:-pipeline}" \
    "$*"
}

brainomics_die() {
  printf '[%s] [%s] error: %s\n' \
    "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    "${BRAINOMICS_STAGE:-pipeline}" \
    "$*" >&2
  exit 1
}

brainomics_is_true() {
  [[ "${1:-}" =~ ^([Tt]|[Tt][Rr][Uu][Ee]|1)$ ]]
}

brainomics_require_command() {
  command -v "$1" >/dev/null 2>&1 || brainomics_die "$1 is unavailable"
}

brainomics_require_file() {
  [ -f "$1" ] || brainomics_die "missing file: $1"
}

brainomics_require_dir() {
  [ -d "$1" ] || brainomics_die "missing directory: $1"
}

brainomics_require_executor() {
  case "$BRAINOMICS_EXECUTOR" in
    local|slurm) ;;
    *) brainomics_die "BRAINOMICS_EXECUTOR must be local or slurm" ;;
  esac
  if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
    brainomics_require_command sbatch
    [ -n "${BRAINOMICS_HPC_PARTITION:-}" ] || brainomics_die "set BRAINOMICS_HPC_PARTITION in hpc/local.env"
  fi
}

# Submit one hpc/*.sbatch script, wait for its exit status, and return its job id.
brainomics_submit() {
  local script="$1"
  shift
  local assignment
  for assignment in "$@"; do
    [[ "$assignment" =~ ^[A-Za-z_][A-Za-z0-9_]*= ]] ||
      brainomics_die "invalid Slurm environment assignment: $assignment"
  done
  local partition="$BRAINOMICS_HPC_PARTITION"
  if grep -q '^# brainomics-partition: gpu' "$script" 2>/dev/null; then
    partition="${BRAINOMICS_HPC_GPU_PARTITION:-$partition}"
  fi
  local log_dir="${BRAINOMICS_HPC_LOG_DIR:-$BRAINOMICS_REPO_ROOT/hpc/logs}"
  mkdir -p "$log_dir"
  local base
  base="$(basename "$script" .sbatch)"
  local sbatch_args=(--parsable --wait --partition="$partition"
    --output="$log_dir/${base}_%j.out" --error="$log_dir/${base}_%j.err")
  if [ -n "${BRAINOMICS_SLURM_MEMORY:-}" ]; then
    sbatch_args+=(--mem="$BRAINOMICS_SLURM_MEMORY")
  fi
  # Set overrides in sbatch's environment, then export it intact. This also
  # preserves spaces/commas in values without Slurm's assignment-list parsing.
  env "$@" sbatch "${sbatch_args[@]}" --export=ALL "$script"
}

brainomics_run_sbatch() {
  local script="$1"
  shift
  brainomics_require_file "$BRAINOMICS_REPO_ROOT/$script"
  local job_id
  brainomics_log "submitting $script and waiting for completion"
  # The submitted job runs the stage body; it must not resubmit itself.
  job_id="$(brainomics_submit "$BRAINOMICS_REPO_ROOT/$script" \
    "BRAINOMICS_EXECUTOR=local" "$@")" ||
    brainomics_die "$script submission or execution failed"
  brainomics_log "$script completed (job $job_id)"
}

# local: run the command here. slurm: submit the matching hpc script instead.
brainomics_dispatch() {
  local sbatch_script="$1"
  shift
  if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
    brainomics_run_sbatch "$sbatch_script" "$@"
  else
    brainomics_log "running locally: $*"
    "$@"
  fi
}
