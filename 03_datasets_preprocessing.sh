#!/usr/bin/env bash
# Reconstruct and standardize each dataset bundle, then refresh the combined
# processed metadata.


set -euo pipefail

BRAINOMICS_STAGE=03_source_preprocessing
source "functions/pipeline_lib.sh"

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
overwrite="${1:-false}"
selected_dataset="${2:-}"

if [ -x "$repo_dir/.venv/bin/python3" ] || [ -x "$repo_dir/.venv/bin/python" ]; then
  export PATH="$repo_dir/.venv/bin:$PATH"
fi

datasets=(
  AllenM1
  EGAS00001006537
  GSE104276
  GSE168408
  GSE186538
  GSE204683
  GSE207334
  GSE212606
  GSE217511
  GSE296073
  GSE67835
  GSE81475
  GSE97942
  HYPOMAP
  Li_et_al_2018
  Ma_et_al_2022
  PRJCA015229
  ROSMAP
  SomaMut
  GSE294786
  Velmeshev_2023
  Wang_2025
)

cd "$repo_dir"

if [ -n "$selected_dataset" ]; then
  supported=false
  for dataset in "${datasets[@]}"; do
    if [ "$dataset" = "$selected_dataset" ]; then
      supported=true
      break
    fi
  done
  if [ "$supported" = false ]; then
    echo "Unsupported formal dataset: $selected_dataset" >&2
    exit 2
  fi
  datasets=("$selected_dataset")
fi

for dataset in "${datasets[@]}"; do
  printf '[%s] reconstructing and standardizing %s\n' \
    "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    "$dataset"
  bash hpc/run_dataset.sh "$dataset" "$overwrite"
done

if [ -z "$selected_dataset" ]; then
  Rscript \
    processing/standardize_datasets.R \
    --combine-only
fi

printf 'Validated %s formal dataset bundle(s)\n' "${#datasets[@]}"
