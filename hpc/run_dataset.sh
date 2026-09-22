#!/usr/bin/env bash

set -euo pipefail

dataset="${1:?Usage: run_dataset.sh <dataset> [overwrite]}"
overwrite="${2:-F}"

case "$dataset" in
  AllenM1|EGAS00001006537|GSE104276|GSE168408|GSE186538|GSE204683|\
  GSE207334|GSE212606|GSE217511|GSE296073|GSE67835|GSE81475|GSE97942|\
  HYPOMAP|Li_et_al_2018|Ma_et_al_2022|PRJCA015229|ROSMAP|SomaMut|\
  GSE294786|Velmeshev_2023|Wang_2025)
    ;;
  *)
    echo "Unsupported dataset: $dataset" >&2
    exit 2
    ;;
esac

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
data_root="$(cd "$repo_dir/../.." && pwd)/data/BrainOmicsData"
source "$repo_dir/functions/utils.sh"
processed_dir="$data_root/processed/$dataset"
run_dir="$processed_dir/run"
mkdir -p "$run_dir"

log_file="$run_dir/processing.log"
exit_file="$run_dir/pipeline.exit"
pid_file="$run_dir/pipeline.pid"
lock_file="$run_dir/pipeline.lock"
rscript_bin="$(command -v Rscript || true)"

if [ -z "$rscript_bin" ] || [ ! -x "$rscript_bin" ]; then
  echo "Rscript is unavailable" >&2
  exit 1
fi

if [ -x "$repo_dir/.venv/bin/python3" ] || [ -x "$repo_dir/.venv/bin/python" ]; then
  export PATH="$repo_dir/.venv/bin:$PATH"
fi

is_true() {
  [[ "$1" =~ ^([Tt]|[Tt][Rr][Uu][Ee]|1)$ ]]
}

is_h5ad_dataset() {
  case "$dataset" in
    Wang_2025|GSE294786|Velmeshev_2023)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

is_original_source_dataset() {
  case "$dataset" in
    AllenM1|EGAS00001006537|GSE144136|\
    GSE178175|GSE202210)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

execute_dataset_pipeline() {
  if [ "$dataset" = "GSE294786" ]; then
    "$rscript_bin" processing/GSE294786.R "$overwrite"
    return
  fi

  if [ "$dataset" = "GSE168408" ]; then
    "$rscript_bin" processing/GSE168408.R "$overwrite"
    return
  fi

  if is_h5ad_dataset; then
    "$rscript_bin" "processing/${dataset}.R" "$overwrite"
    return
  fi

  if is_original_source_dataset; then
    "$rscript_bin" "processing/${dataset}.R" "$overwrite"
  else
    object_file="$processed_dir/${dataset}_processed.rds"
    if is_true "$overwrite" || [ ! -f "$object_file" ]; then
      "$rscript_bin" "processing/${dataset}.R"
    else
      printf 'Retaining existing complete matrix: %s\n' "$object_file"
    fi
  fi

  standardize_args=(
    processing/standardize_datasets.R
    --dataset "$dataset"
    --overwrite
  )
  "$rscript_bin" "${standardize_args[@]}"
}

run_pipeline() {
  printf '%s\n' "$$" > "$pid_file"
  rm -f "$exit_file"
  cd "$repo_dir"

  set +e
  (
    set -e
    execute_dataset_pipeline
  ) 2>&1 | tee "$log_file"
  status="${PIPESTATUS[0]}"
  set -e

  temporary_exit="${exit_file}.tmp.$$"
  printf '%s\n' "$status" > "$temporary_exit"
  mv "$temporary_exit" "$exit_file"
  exit "$status"
}

if command -v flock >/dev/null 2>&1; then
  exec 9>"$lock_file"
  if ! flock -n 9; then
    echo "A processing run is already active for $dataset" >&2
    exit 3
  fi
fi

run_pipeline
