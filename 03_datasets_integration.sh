#!/usr/bin/env bash
set -euo pipefail

source "functions/utils.sh"

check_command Rscript

overwrite="${1:-F}"
integration_dir="../../data/BrainOmicsData/integration"

should_process() {
  local target_file="$1"
  if [[ "$overwrite" =~ ^([Tt]|[Tt][Rr][Uu][Ee]|1)$ ]]; then
    return 0
  fi
  if [ ! -f "$target_file" ]; then
    return 0
  fi
  return 1
}

run_step() {
  local script="$1"
  local label="$2"
  local target_file="$3"

  if should_process "$target_file"; then
    log_message "Running ${label}..."
    Rscript "$script"
    log_message "${label} completed successfully!" --message-type success
  else
    log_message "${label} already completed: ${target_file}"
  fi
}

run_step \
  "integration/datasets_integration_01.R" \
  "integration step 01: metadata harmonization and raw object list" \
  "${integration_dir}/metadata_filtered.rds"

run_step \
  "integration/datasets_integration_02.R" \
  "integration step 02: object merge and gene filtering" \
  "${integration_dir}/objects_filtered.rds"

run_step \
  "integration/datasets_integration_03.R" \
  "integration step 03: PCA, clustering, RPCA/Harmony integration and LISI input" \
  "${integration_dir}/objects_integrated.rds"

run_step \
  "integration/datasets_annotation.R" \
  "integration step 04: cluster-to-cell-type annotation" \
  "${integration_dir}/objects_celltypes.rds"
