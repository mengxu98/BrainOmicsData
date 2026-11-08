#!/usr/bin/env bash
# Rebuild real count matrices from the original data for every source, apply the
# common gene panel and assemble the common matrix.
#
# BRAINOMICS_RUN_ROOT holds inputs/ and the gene-removal policies.
# Source code is read from this checkout; count inputs must be locally available.
set -euo pipefail

BRAINOMICS_STAGE=04_source_counts_assembly
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

stage_args=("$@")
check_inputs=false
selected_dataset=""
while [ "$#" -gt 0 ]; do
  case "$1" in
    --check-inputs) check_inputs=true; shift ;;
    --dataset)
      [ "$#" -ge 2 ] || brainomics_die "--dataset needs a source name"
      selected_dataset="$2"; shift 2 ;;
    *) brainomics_die "Usage: 04_source_counts_assembly.sh [--check-inputs] [--dataset NAME]" ;;
  esac
done

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_run_cpu_stage 04_source_counts_assembly.sh "${stage_args[@]}"
  exit 0
fi

brainomics_require_command "$BRAINOMICS_PYTHON"
"$BRAINOMICS_PYTHON" - "$BRAINOMICS_RUN_ROOT" "$BRAINOMICS_REPO_ROOT" "$selected_dataset" "$check_inputs" <<'PY'
import csv
import os
from pathlib import Path
import sys

run = Path(sys.argv[1])
repo = Path(sys.argv[2])
registry = Path(os.environ.get('BRAINOMICS_DATASET_LIST', str(repo/'data/source_redistribution_evidence.tsv')))
if registry.suffix == '.txt':
    cohort = registry.read_text().splitlines()
else:
    with registry.open() as handle:
        cohort = [row['Dataset'] for row in csv.DictReader(handle, delimiter='\t')]
if len(cohort) != len(set(cohort)) or len(cohort) != 22:
    raise SystemExit('Expected the fixed registry of 22 reference sources')
if sys.argv[3] and sys.argv[3] not in cohort:
    raise SystemExit('Unknown reference source: '+sys.argv[3])
selected = [sys.argv[3]] if sys.argv[3] else cohort
if sys.argv[4] != 'true':
    for dataset in selected:
        destination = run/'inputs/processed'/dataset
        if destination.exists():
            continue
        source = Path(os.environ['BRAINOMICS_DATA_ROOT'])/'processed'/dataset
        schema = source/'processed_object_format.tsv'
        if not schema.is_file():
            raise SystemExit(f'Preprocess the source before staging: {schema}')
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.symlink_to(source.resolve(), target_is_directory=True)
required = [
    'inputs/plan/selected_gene_panel_34387.tsv',
    'inputs/plan/retained_cells_exploration.tsv.gz',
    'code/source_row_audit_map.tsv.gz',
]
if not sys.argv[3]:
    required.extend(['named_gene_panel_policy.json', 'named_gene_panel_24659.tsv'])
    transfer = ['matrices/NAMED_TRANSFER_VERIFIED.json', 'named_transfer_expected.json']
    if any((run / name).exists() for name in transfer):
        required.extend(transfer)
missing = [str(run / name) for name in required if not (run / name).is_file()]
if not missing:
    for dataset in selected:
        directory = run / 'inputs/processed' / dataset
        schema = directory / 'processed_object_format.tsv'
        if not schema.is_file():
            missing.append(str(schema))
            continue
        with schema.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        if len(rows) != 1:
            raise SystemExit(f'Expected one format record: {schema}')
        row = rows[0]
        key = 'Analysis_Object' if row['Format_Type'] == 'paired_full_analysis_objects' else 'Full_Object'
        fields = [key] if sys.argv[3] else [key, 'Canonical_Metadata']
        for field in fields:
            value = row.get(field, '')
            if not value or not (directory / value).is_file():
                missing.append(str(directory / value) if value else f'{schema}: {field}')
legacy_policy = ['applied_gene_row_exclusions.tsv', 'allowed_count_addition.tsv',
                 'applied_gene_removal_policy.json']
if not all((run/name).is_file() for name in legacy_policy):
    references = Path(os.environ.get('BRAINOMICS_GENE_REFERENCE_DIR',
        str(Path(os.environ['BRAINOMICS_DATA_ROOT'])/'raw/gene_references')))
    reference_files = ['hgnc_complete_set.txt', 'Homo_sapiens.gene_info.gz', 'gene_history.gz']
    reference_files += [f'gencode.v{v}.annotation.gtf.gz' for v in (17, 19, 24, 27, 32, 38, 43, 50)]
    if 'AllenM1' in selected:
        reference_files.append('Allen_GRCh38p2_reference.zip')
    if 'HYPOMAP' in selected:
        reference_files.append('Homo_sapiens.GRCh38.98.gtf.gz')
    missing.extend(str(references / name) for name in reference_files if not (references / name).is_file())
    count_features = Path(os.environ['BRAINOMICS_DATA_ROOT'])/'raw/GSE81475/counts.txt'
    if 'GSE81475' in selected and not count_features.is_file():
        missing.append(str(count_features))
if missing:
    raise SystemExit('Missing frozen assembly inputs; stage 03 does not generate the selection policies:\n' +
                     '\n'.join(missing))
print('Assembly input files are present; no matrices were rebuilt or content checks performed')
PY

if [ "$check_inputs" = true ]; then
  exit 0
fi

brainomics_require_command "$BRAINOMICS_RSCRIPT"
build_args=()
if [ -n "$selected_dataset" ]; then build_args=(--dataset "$selected_dataset"); fi
"$BRAINOMICS_PYTHON" functions/run_gene_removal_build.py "${build_args[@]}"
if [ -n "$selected_dataset" ]; then
  brainomics_log "selected source counts rebuilt; common-matrix assembly was not run"
  exit 0
fi
"$BRAINOMICS_PYTHON" functions/run_named_source_assembly.py
"$BRAINOMICS_RSCRIPT" functions/assemble_common_matrix.R "$BRAINOMICS_RUN_ROOT"

brainomics_log "common-gene reference objects written to $BRAINOMICS_RESULTS_DIR"
