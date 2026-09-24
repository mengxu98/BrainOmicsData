#!/usr/bin/env bash
# Rebuild real count matrices from the original data for every source, apply the
# common gene panel and assemble the common matrix.
#
# BRAINOMICS_RUN_ROOT holds inputs/ and the gene-removal policies.
# Source code is read from this checkout; count inputs must be locally available.
#
# Requirements: 512 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=04_source_counts_assembly
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_die "count reconstruction uses local source inputs; set BRAINOMICS_EXECUTOR=local"
fi

brainomics_require_dir "$BRAINOMICS_RUN_ROOT/inputs"
brainomics_require_file "$BRAINOMICS_RUN_ROOT/named_gene_panel_policy.json"
brainomics_require_command "$BRAINOMICS_RSCRIPT"
brainomics_require_command "$BRAINOMICS_PYTHON"

"$BRAINOMICS_PYTHON" functions/run_gene_removal_build.py
"$BRAINOMICS_PYTHON" functions/run_named_source_assembly.py
"$BRAINOMICS_RSCRIPT" functions/assemble_common_matrix.R "$BRAINOMICS_RUN_ROOT"

brainomics_log "common-gene reference objects written to $BRAINOMICS_RESULTS_DIR"
