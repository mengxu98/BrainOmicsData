#!/usr/bin/env bash
# Source this file to select the same environment for every workflow stage.
brainomics_env_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export BRAINOMICS_ENV_ROOT="${BRAINOMICS_ENV_ROOT:-$brainomics_env_repo/.brainomics-env}"
if [ -x "$BRAINOMICS_ENV_ROOT/bin/Rscript" ] && [ -x "$BRAINOMICS_ENV_ROOT/bin/python" ]; then
  export PATH="$BRAINOMICS_ENV_ROOT/bin:$PATH"
  export BRAINOMICS_RSCRIPT="$BRAINOMICS_ENV_ROOT/bin/Rscript"
  export BRAINOMICS_PYTHON="$BRAINOMICS_ENV_ROOT/bin/python"
  # Allow an explicit scVI command while keeping the shared Python as default.
  export BRAINOMICS_SCVI_PYTHON="${BRAINOMICS_SCVI_PYTHON:-$BRAINOMICS_PYTHON}"
  export R_LIBS_USER="$BRAINOMICS_ENV_ROOT/lib/R/site-library"
  export R_LIBS_SITE="$R_LIBS_USER"
  export R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null
  export PYTHONNOUSERSITE=1
  unset R_LIBS PYTHONPATH PYTHONHOME
fi
export BRAINOMICS_RSCRIPT="${BRAINOMICS_RSCRIPT:-Rscript}"
export BRAINOMICS_PYTHON="${BRAINOMICS_PYTHON:-python3}"
export BRAINOMICS_SCVI_PYTHON="${BRAINOMICS_SCVI_PYTHON:-$BRAINOMICS_PYTHON}"
export RETICULATE_PYTHON="$BRAINOMICS_PYTHON"
unset brainomics_env_repo
