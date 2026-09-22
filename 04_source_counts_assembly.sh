#!/usr/bin/env bash
# Rebuild real count matrices from the original data for every source, apply the
# frozen gene panel, assemble the common matrix and accept it.
#
# Requires BRAINOMICS_RUN_ROOT: the run directory that holds inputs/, code/ and
# the frozen gene-removal policies. Runs on the storage host.
#
# Requirements: 512 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=04_source_counts_assembly
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_die "count reconstruction runs on the storage host; use BRAINOMICS_EXECUTOR=local"
fi

[ -n "${BRAINOMICS_RUN_ROOT:-}" ] ||
  brainomics_die "set BRAINOMICS_RUN_ROOT to the run directory holding inputs/ and code/"
export BRAINOMICS_RUN_ROOT
brainomics_require_dir "$BRAINOMICS_RUN_ROOT/inputs"
brainomics_require_file "$BRAINOMICS_RUN_ROOT/named_gene_panel_policy.json"
check_command Rscript
check_command python3

python3 functions/run_gene_removal_build.py
python3 functions/run_named_source_assembly.py
python3 functions/run_final_source_assembly.py
Rscript functions/assemble_common_matrix.R "$BRAINOMICS_RUN_ROOT"

brainomics_log "24,659-gene matrices accepted at $BRAINOMICS_RUN_ROOT/matrices"
