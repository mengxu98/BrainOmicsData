#!/usr/bin/env bash
# Build the per-cell mitochondrial fraction lookup for the released cohort
# from the original source count matrices (see the percent_mito_* scripts here).
#
# The lookup is computed from the original source counts, never from the released
# expression shards. Sources whose feature universe has no MT- gene stay blank
# (percent_mito is written as an empty field, never zero).
#
# Usage: processing/percent_mito_build.sh DATA_ROOT WORK_DIR ANALYSIS_DIR [COHORT_TSV]
#   DATA_ROOT    dataset root holding raw/ and processed/ (default: BRAINOMICS_DATA_ROOT)
#   WORK_DIR     output directory for the lookup and the per-source QC tables
#   ANALYSIS_DIR analysis directory of the frozen run (metadata_working.rds)
#   COHORT_TSV   optional table with Dataset[,Cells]; defaults to the deposit
#                manifest, then to the repository source table
set -euo pipefail

stage_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$stage_dir/.." && pwd)"
source "$repo_dir/functions/log_message.sh"

if [ -x "$repo_dir/.venv/bin/python3" ] || [ -x "$repo_dir/.venv/bin/python" ]; then
  export PATH="$repo_dir/.venv/bin:$PATH"
fi

data_root="${1:-${BRAINOMICS_DATA_ROOT:?set BRAINOMICS_DATA_ROOT}}"
work="${2:?Usage: build.sh DATA_ROOT WORK_DIR ANALYSIS_DIR [COHORT_TSV]}"
analysis="${3:?Usage: build.sh DATA_ROOT WORK_DIR ANALYSIS_DIR [COHORT_TSV]}"
package="${BRAINOMICS_PACKAGE_DIR:-$data_root/ScienceDB}"
cohort="${4:-}"
if [ -z "$cohort" ]; then
  if [ -f "$package/provenance/dataset_manifest.tsv" ]; then
    cohort="$package/provenance/dataset_manifest.tsv"
  else
    cohort="$repo_dir/data/source_access_summary.tsv"
  fi
fi

raw="$data_root/raw"
qc="$work/source_qc"
mkdir -p "$qc"

r_object_datasets=(GSE204683 HYPOMAP Ma_et_al_2022 SomaMut GSE296073 Li_et_al_2018)
matrix_datasets=(GSE217511 PRJCA015229 ROSMAP GSE186538 GSE207334 GSE212606 GSE81475)

if [ -f "$work/percent_mito_by_cell.tsv.gz" ] && [ "${BRAINOMICS_FORCE_PERCENT_MITO:-0}" != "1" ]; then
  log_message "percent_mito lookup already present in $work" --message-type info
  exit 0
fi

log_message "inventorying source feature universes" --message-type info
python3 "$stage_dir/percent_mito_inventory_sources.py" "$raw" "$qc/feature_inventory.tsv"

for dataset in "${r_object_datasets[@]}"; do
  [ -s "$qc/${dataset}_source_qc.tsv.gz" ] && continue
  log_message "source counts: $dataset (R object)" --message-type info
  Rscript "$stage_dir/percent_mito_inventory_r_objects.R" "$raw" "$qc" "$dataset"
done

for dataset in "${matrix_datasets[@]}"; do
  [ -s "$qc/${dataset}_source_qc.tsv.gz" ] && continue
  log_message "source counts: $dataset (matrix)" --message-type info
  Rscript "$stage_dir/percent_mito_calculate_matrix_sources.R" "$raw" "$qc" "$dataset"
done

if [ ! -s "$qc/Wang_2025_source_qc.tsv.gz" ] || [ ! -s "$qc/Velmeshev_2023_source_qc.tsv.gz" ]; then
  log_message "source counts: Wang_2025, Velmeshev_2023 (H5AD raw/X)" --message-type info
  python3 "$stage_dir/percent_mito_calculate_h5ad.py" "$raw" "$qc"
fi

if [ ! -s "$qc/GSE168408_source_qc.tsv.gz" ]; then
  log_message "source counts: GSE168408 (original 10x H5)" --message-type info
  python3 "$stage_dir/percent_mito_calculate_gse168408.py" "$data_root" "$qc"
fi

if [ ! -s "$qc/gse296073_10x/GSE296073_source_qc.tsv.gz" ]; then
  log_message "source counts: GSE296073 (28 original 10x libraries)" --message-type info
  mkdir -p "$qc/gse296073_10x"
  Rscript "$stage_dir/percent_mito_calculate_gse296073.R" "$raw" "$qc/gse296073_10x" GSE296073
fi

log_message "summarizing source inventories" --message-type info
python3 "$stage_dir/percent_mito_summarize_inventory.py" "$qc" "$cohort" "$work/percent_mito_source_inventory.tsv"

log_message "mapping source cells onto the released cohort" --message-type info
Rscript "$stage_dir/percent_mito_map_to_cohort.R" "$analysis" "$data_root" "$qc" "$work"

python3 - "$work/mapping_status.json" <<'PY'
import json, sys
state = json.load(open(sys.argv[1]))
assert state['state'] == 'LOOKUP_READY', state['state']
print('percent_mito lookup ready:', state['cells'], 'cells,', state['nonmissing'], 'with a source value')
PY
