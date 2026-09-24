#!/usr/bin/env bash
# Verify the single R/Python environment; --install restores it first.
#
# Requires the pinned R and Python environments in environment/ and network
# access when they are installed rather than verified.
set -euo pipefail

BRAINOMICS_STAGE=01_environment
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

case "${1:---verify}" in
  --verify) ;;
  --install)
    if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
      brainomics_run_sbatch hpc/install_environment.sbatch
    else
      bash environment/restore.sh
    fi
    ;;
  *) echo 'Usage: 01_environment.sh [--verify|--install]' >&2; exit 2 ;;
esac
source environment/activate.sh
check_command "$BRAINOMICS_RSCRIPT"
check_command "$BRAINOMICS_PYTHON"
for command in pigz pdfcrop pdfinfo xelatex fc-match; do
  check_command "$command"
done
font_match="$(fc-match Arial --format '%{family}\n' | head -n 1)"
if ! printf '%s\n' "$font_match" | tr ',' '\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | grep -Fxq Arial; then
  brainomics_die "Arial is unavailable to fontconfig (matched: ${font_match:-none})"
fi
"$BRAINOMICS_RSCRIPT" --vanilla environment/verify_environment.R "$BRAINOMICS_REPO_ROOT"
"$BRAINOMICS_PYTHON" environment/verify_python.py --imports

brainomics_log "Unified R 4.5.1 and Python environment match the locked versions"
