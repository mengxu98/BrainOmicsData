#!/bin/bash
#SBATCH -J brainomics_integration_deps
#SBATCH -p cpu-partition
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/install_integration_deps_%j.out
#SBATCH -e /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/install_integration_deps_%j.err

set -euo pipefail

unset LC_ALL LC_CTYPE LANGUAGE CPLUS_INCLUDE_PATH C_INCLUDE_PATH
export LANG=C
module load mathlib/hdf5/1.8.20-gcc

repo_dir="/path/to/hpc/home/repositories/BrainOmicsData"
r_prefix="/path/to/hpc/home/miniforge3"
rscript="$r_prefix/bin/Rscript"
r_libs_user="${R_LIBS_USER:-/path/to/hpc/home/R/library/4.5}"
source_dir="/path/to/hpc/home/software_sources/r-cran"

export R_LIBS_USER="$r_libs_user"
export PATH="$r_prefix/bin:$PATH"
export PKG_CONFIG_PATH="$r_prefix/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export CPATH="$r_prefix/include:${CPATH:-}"
export LIBRARY_PATH="$r_prefix/lib:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$r_prefix/lib:${LD_LIBRARY_PATH:-}"

mkdir -p "$repo_dir/hpc/logs" "$R_LIBS_USER"

fetch_pinned_source_package() {
  local source_url="$1"
  local tarball_name="$2"
  local tarball_sha256="$3"
  local tarball="$source_dir/$tarball_name"
  local temporary="${tarball}.tmp.$$"

  mkdir -p "$source_dir"
  if [ ! -f "$tarball" ]; then
    curl \
      --fail \
      --location \
      --retry 3 \
      --output "$temporary" \
      "$source_url"
    printf '%s  %s\n' "$tarball_sha256" "$temporary" |
      sha256sum -c -
    mv "$temporary" "$tarball"
  fi
  printf '%s  %s\n' "$tarball_sha256" "$tarball" |
    sha256sum -c -
}

install_pinned_source_package() {
  local package_name="$1"
  local package_version="$2"
  local tarball_name="$3"
  local tarball_sha256="$4"
  local tarball="$source_dir/$tarball_name"

  if "$rscript" -e "quit(status = if (requireNamespace('$package_name', quietly = TRUE) && packageVersion('$package_name') == package_version('$package_version')) 0 else 1)"; then
    printf '%s %s already installed\n' "$package_name" "$package_version"
    return 0
  fi
  if [ ! -f "$tarball" ]; then
    echo "Missing pinned R source tarball: $tarball" >&2
    exit 1
  fi
  printf '%s  %s\n' "$tarball_sha256" "$tarball" |
    sha256sum -c -
  "$r_prefix/bin/R" CMD INSTALL \
    --library="$R_LIBS_USER" \
    "$tarball"
}

install_pinned_lisi() {
  local package_version="1.0"
  local tarball_name="lisi_a917556310d8d2c66833dcc35aa3d0f4d1b6e0f4.tar.gz"
  local tarball_sha256="b4724c19c1bbcc7891f63ce70cae61da0b8870819981a8a40e864ed1e3caf884"
  local tarball="$source_dir/$tarball_name"
  local numeric_patch="$repo_dir/hpc/patches/lisi_double_precision.patch"
  local numeric_patch_sha256="6c8d1e25cf62aaa385ef9fff90eeedfb31519768a9853c5d5a116a85837171dd"
  local marker_file="$R_LIBS_USER/lisi/BRAINOMICS_PINNED_SOURCE"
  local stage_root package_dir makevars_file marker_temporary

  if "$rscript" -e "quit(status = if (requireNamespace('lisi', quietly = TRUE) && packageVersion('lisi') == package_version('$package_version')) 0 else 1)" &&
     [ -f "$marker_file" ] &&
     grep -Fxq "source_commit=a917556310d8d2c66833dcc35aa3d0f4d1b6e0f4" "$marker_file" &&
     grep -Fxq "tarball_sha256=$tarball_sha256" "$marker_file" &&
     grep -Fxq "build_patch=CXX11_to_CXX14" "$marker_file" &&
     grep -Fxq "numeric_patch=double_precision_probability_normalization_v1" "$marker_file" &&
     grep -Fxq "numeric_patch_sha256=$numeric_patch_sha256" "$marker_file"; then
    printf 'lisi %s already installed\n' "$package_version"
    return 0
  fi
  printf '%s  %s\n' "$tarball_sha256" "$tarball" |
    sha256sum -c -
  stage_root="$(mktemp -d "${TMPDIR:-/tmp}/brainomics-lisi.XXXXXX")"
  tar -xzf "$tarball" -C "$stage_root"
  package_dir="$(
    find "$stage_root" \
      -mindepth 1 \
      -maxdepth 1 \
      -type d \
      -name 'LISI-*' \
      -print \
      -quit
  )"
  if [ -z "$package_dir" ]; then
    echo "Pinned LISI source archive has an unexpected layout" >&2
    exit 1
  fi
  makevars_file="$package_dir/src/Makevars"
  "$rscript" - "$makevars_file" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
file <- args[[1L]]
lines <- readLines(file, warn = FALSE)
target <- "CXX_STD = CXX11"
if (sum(lines == target) != 1L) {
  stop("Pinned LISI Makevars no longer has the reviewed CXX11 declaration")
}
  lines[lines == target] <- "CXX_STD = CXX14"
  writeLines(lines, file, useBytes = TRUE)
RS
  if [ ! -f "$numeric_patch" ]; then
    echo "Missing reviewed LISI numerical patch: $numeric_patch" >&2
    exit 1
  fi
  printf '%s  %s\n' "$numeric_patch_sha256" "$numeric_patch" |
    sha256sum -c -
  patch --batch --forward -d "$package_dir" -p1 < "$numeric_patch"
  # LISI 1.0 predates current RcppArmadillo, which now requires C++14. The
  # build compatibility change and reviewed double-precision numerical patch
  # are applied to the verified temporary source tree. The latter prevents
  # single-precision underflow in probability normalization for isolated cells.
  "$r_prefix/bin/R" CMD INSTALL \
    --library="$R_LIBS_USER" \
    "$package_dir"
  marker_temporary="${marker_file}.tmp.$$"
  {
    printf '%s\n' \
      "source_commit=a917556310d8d2c66833dcc35aa3d0f4d1b6e0f4" \
      "tarball_sha256=$tarball_sha256" \
      "build_patch=CXX11_to_CXX14" \
      "numeric_patch=double_precision_probability_normalization_v1" \
      "numeric_patch_sha256=$numeric_patch_sha256"
  } > "$marker_temporary"
  mv "$marker_temporary" "$marker_file"
  rm -rf -- "$stage_root"
}

install_pinned_source_package \
  RhpcBLASctl 0.23-42 RhpcBLASctl_0.23-42.tar.gz \
  5c889d5b69e264060b9f1f0383c447f594855b8afc15b7d76d39e4d62b946615
fetch_pinned_source_package \
  "https://cran.r-project.org/src/contrib/R.methodsS3_1.8.2.tar.gz" \
  "R.methodsS3_1.8.2.tar.gz" \
  "822d5e61dad4c91e8883be2b38d7b89f87492046d0fe345704eb5d2658927c2e"
install_pinned_source_package \
  R.methodsS3 1.8.2 R.methodsS3_1.8.2.tar.gz \
  822d5e61dad4c91e8883be2b38d7b89f87492046d0fe345704eb5d2658927c2e
fetch_pinned_source_package \
  "https://cran.r-project.org/src/contrib/R.oo_1.27.1.tar.gz" \
  "R.oo_1.27.1.tar.gz" \
  "5faf599cc1f027d8b80e1270aa6a43e29dccd6fee8287953f13e7ab3d94270f5"
install_pinned_source_package \
  R.oo 1.27.1 R.oo_1.27.1.tar.gz \
  5faf599cc1f027d8b80e1270aa6a43e29dccd6fee8287953f13e7ab3d94270f5
fetch_pinned_source_package \
  "https://cran.r-project.org/src/contrib/R.utils_2.13.0.tar.gz" \
  "R.utils_2.13.0.tar.gz" \
  "ab2043c34e129928ff85a037ce7a3f7791f245f49a55ee9a95bd7df0813bcf25"
install_pinned_source_package \
  R.utils 2.13.0 R.utils_2.13.0.tar.gz \
  ab2043c34e129928ff85a037ce7a3f7791f245f49a55ee9a95bd7df0813bcf25
install_pinned_source_package \
  harmony 2.0.5 harmony_2.0.5.tar.gz \
  907a3c4808656f6bad5a8e314704f92a0cb159ae298390315d5c73794d9a0c6c
fetch_pinned_source_package \
  "https://codeload.github.com/immunogenomics/LISI/tar.gz/a917556310d8d2c66833dcc35aa3d0f4d1b6e0f4" \
  "lisi_a917556310d8d2c66833dcc35aa3d0f4d1b6e0f4.tar.gz" \
  "b4724c19c1bbcc7891f63ce70cae61da0b8870819981a8a40e864ed1e3caf884"
install_pinned_lisi

"$rscript" - <<'RS'
stopifnot(
  requireNamespace("RhpcBLASctl", quietly = TRUE),
  packageVersion("RhpcBLASctl") == package_version("0.23-42"),
  requireNamespace("R.methodsS3", quietly = TRUE),
  packageVersion("R.methodsS3") == package_version("1.8.2"),
  requireNamespace("R.oo", quietly = TRUE),
  packageVersion("R.oo") == package_version("1.27.1"),
  requireNamespace("R.utils", quietly = TRUE),
  packageVersion("R.utils") == package_version("2.13.0"),
  requireNamespace("harmony", quietly = TRUE),
  packageVersion("harmony") == package_version("2.0.5"),
  requireNamespace("lisi", quietly = TRUE),
  packageVersion("lisi") == package_version("1.0"),
  requireNamespace("Seurat", quietly = TRUE),
  "HarmonyIntegration" %in% getNamespaceExports("Seurat")
)
lisi_pin <- readLines(file.path(
  find.package("lisi"),
  "BRAINOMICS_PINNED_SOURCE"
))
stopifnot(identical(
  lisi_pin,
  c(
    "source_commit=a917556310d8d2c66833dcc35aa3d0f4d1b6e0f4",
    "tarball_sha256=b4724c19c1bbcc7891f63ce70cae61da0b8870819981a8a40e864ed1e3caf884",
    "build_patch=CXX11_to_CXX14",
    "numeric_patch=double_precision_probability_normalization_v1",
    "numeric_patch_sha256=6c8d1e25cf62aaa385ef9fff90eeedfb31519768a9853c5d5a116a85837171dd"
  )
))

set.seed(20260730L)
make_fixture <- function(project) {
  counts <- Matrix::rsparsematrix(
    nrow = 100L,
    ncol = 24L,
    density = 0.25,
    rand.x = function(n) stats::rpois(n, lambda = 2) + 1
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0(project, "_cell", seq_len(ncol(counts)))
  SeuratObject::CreateSeuratObject(
    counts = counts,
    project = project
  )
}
fixture <- merge(
  make_fixture("batch_a"),
  make_fixture("batch_b")
)
fixture <- Seurat::NormalizeData(fixture, verbose = FALSE)
fixture <- Seurat::FindVariableFeatures(
  fixture,
  nfeatures = 30L,
  verbose = FALSE
)
fixture <- Seurat::ScaleData(
  fixture,
  features = SeuratObject::VariableFeatures(fixture),
  verbose = FALSE
)
fixture <- Seurat::RunPCA(
  fixture,
  features = SeuratObject::VariableFeatures(fixture),
  npcs = 10L,
  seed.use = 20260730L,
  verbose = FALSE
)
fixture <- Seurat::IntegrateLayers(
  object = fixture,
  method = Seurat::HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  verbose = FALSE
)
stopifnot(
  "integrated.harmony" %in% names(fixture@reductions),
  nrow(Seurat::Embeddings(fixture, "integrated.harmony")) ==
    ncol(fixture)
)
lisi_fixture <- lisi::compute_lisi(
  X = Seurat::Embeddings(fixture, "pca"),
  meta_data = data.frame(
    Dataset = fixture$orig.ident,
    row.names = colnames(fixture)
  ),
  label_colnames = "Dataset",
  perplexity = 10
)
stopifnot(
  identical(names(lisi_fixture), "Dataset"),
  nrow(lisi_fixture) == ncol(fixture),
  all(is.finite(lisi_fixture$Dataset)),
  all(lisi_fixture$Dataset >= 1 - 1e-12)
)
extreme_distances <- matrix(seq(104, 109, length.out = 89L), ncol = 1L)
extreme_indices <- matrix(as.integer(0:88), ncol = 1L)
extreme_labels <- rep(0, 89L)
extreme_lisi <- 1 / lisi::compute_simpson_index(
  extreme_distances,
  extreme_indices,
  extreme_labels,
  n_batches = 1L,
  perplexity = 30,
  tol = 1e-5
)
stopifnot(
  length(extreme_lisi) == 1L,
  is.finite(extreme_lisi),
  abs(extreme_lisi - 1) < 1e-12
)
cat(
  "RhpcBLASctl",
  as.character(packageVersion("RhpcBLASctl")),
  "and harmony",
  as.character(packageVersion("harmony")),
  "and lisi",
  as.character(packageVersion("lisi")),
  "validated\n"
)
RS
