#!/usr/bin/env bash
# Shared paths, logging and HPC submission helpers for the numbered
# drivers 01_-12_ in the repository root.

set -euo pipefail

BRAINOMICS_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BRAINOMICS_REPO_ROOT="$(cd "$BRAINOMICS_LIB_DIR/.." && pwd)"
BRAINOMICS_DATA_ROOT="${BRAINOMICS_DATA_ROOT:-$(cd "$BRAINOMICS_REPO_ROOT/../.." && pwd)/data/BrainOmicsData}"
BRAINOMICS_RESULTS_DIR="${BRAINOMICS_RESULTS_DIR:-$BRAINOMICS_DATA_ROOT/integration_25}"
BRAINOMICS_EXECUTOR="${BRAINOMICS_EXECUTOR:-local}"
BRAINOMICS_POLL_SECONDS="${BRAINOMICS_POLL_SECONDS:-60}"

export BRAINOMICS_REPO_ROOT BRAINOMICS_DATA_ROOT BRAINOMICS_RESULTS_DIR
export BRAINOMICS_EXECUTOR

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
    local|hpc) ;;
    *) brainomics_die "BRAINOMICS_EXECUTOR must be local or hpc" ;;
  esac
  if [ "$BRAINOMICS_EXECUTOR" = "hpc" ]; then
    brainomics_require_command sbatch
    brainomics_require_command squeue
  fi
}

# Submit one hpc/*.sbatch script and return its job id in stdout.
brainomics_submit() {
  local script="$1"
  shift
  local export_spec="ALL"
  if [ "$#" -gt 0 ]; then
    export_spec="ALL,$*"
  fi
  sbatch --parsable --export="$export_spec" "$script"
}

# Block until the job leaves the queue, then report its accounting state.
brainomics_wait_job() {
  local job_id="$1"
  local state=""
  while squeue -h -j "$job_id" | grep -q .; do
    sleep "$BRAINOMICS_POLL_SECONDS"
  done
  if command -v sacct >/dev/null 2>&1; then
    state="$(sacct -n -X -j "$job_id" --format=State 2>/dev/null | head -n 1 | tr -d ' ')"
  fi
  case "$state" in
    ""|COMPLETED|COMPLETING)
      return 0
      ;;
    *)
      brainomics_die "job $job_id finished with state $state"
      ;;
  esac
}

brainomics_run_sbatch() {
  local script="$1"
  shift
  brainomics_require_file "$BRAINOMICS_REPO_ROOT/$script"
  local job_id
  # The submitted job runs the stage body; it must not resubmit itself.
  job_id="$(brainomics_submit "$BRAINOMICS_REPO_ROOT/$script" \
    "BRAINOMICS_EXECUTOR=local" "$@")"
  brainomics_log "submitted $script as job $job_id"
  brainomics_wait_job "$job_id"
  brainomics_log "$script completed (job $job_id)"
}

# local: run the command here. hpc: submit the matching hpc script instead.
brainomics_dispatch() {
  local sbatch_script="$1"
  shift
  if [ "$BRAINOMICS_EXECUTOR" = "hpc" ]; then
    brainomics_run_sbatch "$sbatch_script" "$@"
  else
    brainomics_log "running locally: $*"
    "$@"
  fi
}
