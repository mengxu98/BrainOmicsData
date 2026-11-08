#!/usr/bin/env bash
# Restore the unified Linux environment from public package repositories.
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[ "$#" -eq 0 ] || { echo 'Usage: restore.sh' >&2; exit 2; }
if [ "$(uname -s)" != Linux ] || [ "$(uname -m)" != x86_64 ]; then
  echo 'The canonical CUDA environment targets Linux x86_64.' >&2
  exit 2
fi
export BRAINOMICS_ENV_ROOT="${BRAINOMICS_ENV_ROOT:-$repo_dir/.brainomics-env}"
micromamba_bin="${BRAINOMICS_MICROMAMBA:-micromamba}"
command -v "$micromamba_bin" >/dev/null 2>&1 || {
  echo 'Install micromamba or set BRAINOMICS_MICROMAMBA to its executable; see environment/REPRODUCIBILITY.md.' >&2
  exit 1
}
# One shared prefix, including compilers, rendering tools and libraries needed
# by the workflow.
platform_packages=(
  'r-base=4.5.1' 'python=3.12.8' pip make compilers cmake pkg-config
  curl git patch hdf5 'icu=75.1' 'libxml2=2.15.1' 'libxml2-devel=2.15.1'
  libcurl openssl expat fontconfig freetype
  harfbuzz fribidi libpng libtiff libjpeg-turbo cairo zlib pandoc poppler
  pigz ghostscript texlive-core xorg-xorgproto liblzma-devel
)
if [ ! -f "$BRAINOMICS_ENV_ROOT/conda-meta/history" ]; then
  "$micromamba_bin" create --yes --prefix "$BRAINOMICS_ENV_ROOT" \
    --override-channels --channel conda-forge --strict-channel-priority \
    "${platform_packages[@]}"
else
  "$micromamba_bin" install --yes --prefix "$BRAINOMICS_ENV_ROOT" \
    --override-channels --channel conda-forge --strict-channel-priority \
    "${platform_packages[@]}"
fi
export R_LIBS_USER="$BRAINOMICS_ENV_ROOT/lib/R/site-library"
mkdir -p "$R_LIBS_USER" "$BRAINOMICS_ENV_ROOT/brainomics"
export R_MAKEVARS_USER="$BRAINOMICS_ENV_ROOT/brainomics/Makevars"
touch "$R_MAKEVARS_USER"
export BRAINOMICS_RSCRIPT="$BRAINOMICS_ENV_ROOT/bin/Rscript"
export BRAINOMICS_PYTHON="$BRAINOMICS_ENV_ROOT/bin/python"
export BRAINOMICS_SCVI_PYTHON="$BRAINOMICS_PYTHON"
source "$repo_dir/environment/activate.sh"
export HOST=x86_64-conda-linux-gnu
export CC="$BRAINOMICS_ENV_ROOT/bin/$HOST-cc"
export PKG_CONFIG_LIBDIR="$BRAINOMICS_ENV_ROOT/lib/pkgconfig:$BRAINOMICS_ENV_ROOT/share/pkgconfig"
export PKG_CONFIG_PATH="$PKG_CONFIG_LIBDIR"
pkg-config --cflags --libs fontconfig freetype2 libxml-2.0 hdf5 >/dev/null

# These checks precede package installation even when reusing an existing prefix.
"$micromamba_bin" run --prefix "$BRAINOMICS_ENV_ROOT" \
  Rscript --vanilla -e 'stopifnot(identical(as.character(getRversion()), "4.5.1"))'
"$micromamba_bin" run --prefix "$BRAINOMICS_ENV_ROOT" \
  python -c 'import sys; assert sys.version_info[:3] == (3, 12, 8), sys.version'

"$micromamba_bin" run --prefix "$BRAINOMICS_ENV_ROOT" \
  Rscript --vanilla "$repo_dir/environment/restore_r.R" "$repo_dir"
"$micromamba_bin" run --prefix "$BRAINOMICS_ENV_ROOT" \
  Rscript --vanilla "$repo_dir/environment/verify_environment.R" "$repo_dir"

if ! command -v xelatex >/dev/null 2>&1 || ! kpsewhich fontspec.sty >/dev/null 2>&1 ||
   ! kpsewhich tikz.sty >/dev/null 2>&1 || ! kpsewhich geometry.sty >/dev/null 2>&1; then
  "$BRAINOMICS_RSCRIPT" --vanilla - <<'RS'
target <- file.path(Sys.getenv("BRAINOMICS_ENV_ROOT"), "brainomics", "TinyTeX")
Sys.setenv(TEXLIVE_DOWNLOADER = "curl")
archive <- file.path(dirname(target), "TinyTeX-1-v2023.04.tar.gz")
if (!file.exists(archive)) {
  partial <- paste0(archive, ".part")
  stopifnot(system2("curl", c("--fail", "--location", "--silent", "--show-error",
    "--retry", "3", "--retry-all-errors", "--output", shQuote(partial),
    "https://github.com/rstudio/tinytex-releases/releases/download/v2023.04/TinyTeX-1-v2023.04.tar.gz")) == 0L)
  stopifnot(file.rename(partial, archive))
}
tinytex::install_tinytex(force = TRUE, dir = target, version = "2023.04",
                         bundle = archive, repository = "ctan", extra_packages = character(), add_path = FALSE)
binary <- Sys.glob(file.path(target, "bin", "*", "tlmgr"))
stopifnot(length(binary) == 1L)
stopifnot(system2(binary, c("option", "repository", "https://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2023/tlnet-final")) == 0L)
status <- system2(binary, c("update", "--self"))
stopifnot(status %in% c(0L, 100L))
stopifnot(system2(binary, c("install", "xetex", "latex-bin", "fontspec", "geometry", "pgf", "graphics", "pdfcrop")) == 0L)
RS
  source "$repo_dir/environment/activate.sh"
fi
for tex_input in article.cls fontspec.sty geometry.sty graphicx.sty tikz.sty; do
  kpsewhich "$tex_input" >/dev/null || { echo "Missing TeX input: $tex_input" >&2; exit 1; }
done

# PyTorch's CUDA 12.4 wheel is selected explicitly, independently of the host
# driver's version. A working NVIDIA driver is required for GPU training.
"$micromamba_bin" run --prefix "$BRAINOMICS_ENV_ROOT" python -m pip install \
  --index-url "${PIP_INDEX_URL:-https://pypi.org/simple}" \
  --extra-index-url https://download.pytorch.org/whl/cu124 \
  --requirement "$repo_dir/environment/python-scvi-freeze.txt"
"$micromamba_bin" run --prefix "$BRAINOMICS_ENV_ROOT" python -m pip check
"$micromamba_bin" run --prefix "$BRAINOMICS_ENV_ROOT" \
  python "$repo_dir/environment/verify_python.py" --imports
"$micromamba_bin" list --prefix "$BRAINOMICS_ENV_ROOT" --explicit \
  > "$BRAINOMICS_ENV_ROOT/brainomics/platform-explicit.txt"
echo "Environment restored at $BRAINOMICS_ENV_ROOT"
echo "Activate it with: source $repo_dir/environment/activate.sh"
