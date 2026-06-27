#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ENV_NAME:-brainomics-h5ad}"
CONDA_BIN="${CONDA_BIN:-/work/home/mengxu1310/miniforge3/bin/conda}"
MAMBA_BIN="${MAMBA_BIN:-/work/home/mengxu1310/miniforge3/bin/mamba}"
SCOP_SOURCE_DIR="${SCOP_SOURCE_DIR:-/work/home/mengxu1310/repositories/scop}"
CRAN_REPO="${CRAN_REPO:-https://mirrors.tuna.tsinghua.edu.cn/CRAN}"
REPO_DIR="${REPO_DIR:-/work/home/mengxu1310/repositories/BrainOmicsData}"
INSTALL_FULL_SCOP="${INSTALL_FULL_SCOP:-0}"

ENV_PREFIX="$("$CONDA_BIN" env list | awk -v env="$ENV_NAME" '$1 == env {print $NF}')"
if [[ -z "$ENV_PREFIX" || ! -d "$ENV_PREFIX" ]]; then
  "$MAMBA_BIN" create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    python=3.11 numpy scipy h5py anndata \
    r-base=4.5 r-seurat r-seuratobject r-reticulate r-matrix \
    r-ggplot2 r-dplyr r-igraph r-uwot r-cli r-rlang r-gtable \
    r-rcpp r-rcppparallel r-ggrepel r-ggforce r-ggnewscale r-pak \
    bioconductor-complexheatmap bioconductor-genomicranges \
    bioconductor-s4vectors bioconductor-summarizedexperiment r-signac
  ENV_PREFIX="$("$CONDA_BIN" env list | awk -v env="$ENV_NAME" '$1 == env {print $NF}')"
else
  echo "Conda environment already exists: $ENV_PREFIX"
fi

if [[ -z "$ENV_PREFIX" || ! -d "$ENV_PREFIX" ]]; then
  echo "Failed to locate conda environment: $ENV_NAME" >&2
  exit 1
fi

export RETICULATE_PYTHON="$ENV_PREFIX/bin/python"

if ! "$ENV_PREFIX/bin/Rscript" -e 'quit(status = !requireNamespace("pak", quietly = TRUE))'; then
  "$MAMBA_BIN" install -y -n "$ENV_NAME" -c conda-forge r-pak
fi

rm -rf "$ENV_PREFIX/lib/R/library"/00LOCK*

if [[ "$INSTALL_FULL_SCOP" == "1" ]]; then
"$ENV_PREFIX/bin/Rscript" - <<RS
options(repos = c(CRAN = "$CRAN_REPO"))
missing <- c("thisplot", "thisutils")[!vapply(c("thisplot", "thisutils"), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)
RS
fi

if [[ "$INSTALL_FULL_SCOP" == "1" && -d "$SCOP_SOURCE_DIR" ]]; then
  "$ENV_PREFIX/bin/R" CMD INSTALL "$SCOP_SOURCE_DIR" || {
    echo "Full scop installation failed; installing lightweight scop h5ad exporter." >&2
    "$ENV_PREFIX/bin/Rscript" "$REPO_DIR/hpc/install_lightweight_scop_h5ad.R"
  }
elif [[ "$INSTALL_FULL_SCOP" == "1" ]]; then
  echo "Missing scop source directory: $SCOP_SOURCE_DIR" >&2
  echo "Copy /Users/mx/Study/repositories/scop to $SCOP_SOURCE_DIR and rerun." >&2
  exit 1
else
  "$ENV_PREFIX/bin/Rscript" "$REPO_DIR/hpc/install_lightweight_scop_h5ad.R"
fi

"$ENV_PREFIX/bin/Rscript" - <<'RS'
stopifnot(requireNamespace("scop", quietly = TRUE))
stopifnot(exists("srt_to_h5ad", envir = asNamespace("scop"), inherits = FALSE))
stopifnot(reticulate::py_module_available("anndata"))
stopifnot(reticulate::py_module_available("numpy"))
cat("brainomics-h5ad environment is ready.\n")
RS
