# Restore one R library for preprocessing, integration, evaluation and plotting.
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args)) args[[1L]] else ".")
lock <- read.delim(file.path(repo, "environment", "r-packages.lock.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
lock$Archive_SHA256[is.na(lock$Archive_SHA256) | lock$Archive_SHA256 == "-"] <- ""
stopifnot(!anyDuplicated(lock$Package), sum(lock$Package == "R") == 1L,
          identical(as.character(getRversion()), lock$Version[lock$Package == "R"]))
lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib)) stop("Set R_LIBS_USER to the shared environment library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(lib), .Library))
options(timeout = max(600L, getOption("timeout")),
        BioC_mirror = "https://bioconductor.posit.co",
        repos = c(CRAN = "https://cran.r-project.org",
                  BioCsoft = "https://bioconductor.posit.co/packages/3.22/bioc",
                  BioCann = "https://bioconductor.posit.co/packages/3.22/data/annotation",
                  BioCexp = "https://bioconductor.posit.co/packages/3.22/data/experiment"))
Sys.setenv(R_BIOC_VERSION = "3.22", PKG_SYSREQS = "false",
           PKG_SYSREQS_DB_UPDATE = "false", ENABLE_INSTALL_COUNTING = "no")
work <- tempfile("brainomics-r-restore-")
dir.create(work)
# This is a script, so a top-level on.exit() is not available. Temporary source
# trees remain outside the repository if installation stops unexpectedly.

sha256 <- function(path) {
  command <- Sys.which("sha256sum")
  if (!nzchar(command)) stop("sha256sum is required to verify source archives")
  output <- system2(command, shQuote(path), stdout = TRUE)
  if (!is.null(attr(output, "status"))) stop("Cannot hash ", path)
  sub(" .*", "", output[[1L]])
}
download <- function(urls, filename, expected = "") {
  destination <- file.path(work, filename)
  success <- FALSE
  for (url in rep(urls, 3L)) {
    success <- tryCatch({
      status <- utils::download.file(url, destination, mode = "wb", quiet = TRUE)
      identical(status, 0L) && file.info(destination)$size > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (success) break
  }
  if (!success) stop("Cannot download pinned source: ", paste(urls, collapse = ", "))
  if (nzchar(expected) && !identical(sha256(destination), expected)) {
    stop("Source archive checksum mismatch: ", filename)
  }
  destination
}
cran_tarball <- function(package, version, expected = "") {
  filename <- paste0(package, "_", version, ".tar.gz")
  download(c(paste0("https://cran.r-project.org/src/contrib/", filename),
             paste0("https://cran.r-project.org/src/contrib/Archive/", package, "/", filename)),
           filename, expected)
}
unpack <- function(archive, folder, package = NULL) {
  target <- file.path(work, folder)
  dir.create(target)
  utils::untar(archive, exdir = target)
  descriptions <- list.files(target, pattern = "^DESCRIPTION$", recursive = TRUE, full.names = TRUE)
  if (!is.null(package)) {
    descriptions <- descriptions[vapply(descriptions, function(path) {
      identical(unname(read.dcf(path)[1L, "Package"]), package)
    }, logical(1L))]
  }
  if (length(descriptions) != 1L) stop("Unexpected source archive layout: ", archive)
  dirname(descriptions[[1L]])
}

# pak is self-contained and bootstraps using only R's standard packages.
pak_version <- lock$Version[lock$Package == "pak"]
if (length(pak_version) != 1L) stop("The installer package pak must have one locked version")
current_pak <- tryCatch(utils::packageVersion("pak"), error = function(e) NULL)
if (is.null(current_pak) || current_pak != base::package_version(pak_version)) {
  utils::install.packages(cran_tarball("pak", pak_version), repos = NULL,
                          lib = lib, type = "source")
}
if (!requireNamespace("pak", quietly = TRUE) ||
    utils::packageVersion("pak") != base::package_version(pak_version)) {
  stop("The locked pak installer could not be restored")
}

# Install the pinned LISI source and numerical patch.
pin_lines <- readLines(file.path(repo, "environment", "lisi", "LISI_PINNED_SOURCE"))
pin <- setNames(sub("^[^=]+=", "", pin_lines), sub("=.*$", "", pin_lines))
lisi_archive <- download(
  paste0("https://codeload.github.com/immunogenomics/LISI/tar.gz/", pin[["source_commit"]]),
  "lisi.tar.gz", pin[["tarball_sha256"]])
lisi_source <- unpack(lisi_archive, "lisi")
makevars <- file.path(lisi_source, "src", "Makevars")
lines <- readLines(makevars)
if (sum(lines == "CXX_STD = CXX11") != 1L) stop("Unexpected LISI Makevars")
lines[lines == "CXX_STD = CXX11"] <- "CXX_STD = CXX14"
writeLines(lines, makevars)
patch <- file.path(repo, "hpc", "patches", "lisi_double_precision.patch")
stopifnot(identical(sha256(patch), pin[["numeric_patch_sha256"]]))
if (system2("patch", c("--batch", "--forward", "-d", shQuote(lisi_source), "-p1"),
            stdin = patch) != 0L) stop("Unable to apply the pinned LISI patch")

scop_pin <- read.delim(file.path(repo, "environment", "scop-plotting.lock.tsv"),
                      stringsAsFactors = FALSE)
scop_archive <- file.path(repo, "environment", scop_pin$Source_Archive)
stopifnot(nrow(scop_pin) == 1L, identical(sha256(scop_archive), scop_pin$Archive_SHA256),
          identical(scop_pin$Version, lock$Version[lock$Package == "scop"]))
scop_source <- unpack(scop_archive, "scop")
desc <- read.dcf(file.path(scop_source, "DESCRIPTION"))
desc <- desc[, colnames(desc) != "Remotes", drop = FALSE]
write.dcf(desc, file.path(scop_source, "DESCRIPTION"))

packages <- lock[lock$Package != "R" & lock$Source != "R base", , drop = FALSE]
installed <- utils::installed.packages(lib.loc = .libPaths())
present <- match(packages$Package, installed[, "Package"])
same_version <- !is.na(present) &
  installed[present, "Version"] == packages$Version
source_pinned <- packages$Source == "GitHub" |
  packages$Package %in% c("lisi", "scop") | nzchar(packages$Archive_SHA256)
packages <- packages[!same_version | source_pinned, , drop = FALSE]
refs <- vapply(seq_len(nrow(packages)), function(i) {
  entry <- packages[i, ]
  name <- entry$Package
  if (name == "lisi") return(paste0("lisi=local::", lisi_source))
  if (name == "scop") return(paste0("scop=local::", scop_source))
  if (entry$Source == "GitHub") {
    reference <- strsplit(entry$Reference, "@", fixed = TRUE)[[1L]]
    stopifnot(length(reference) == 2L)
    archive <- download(paste0("https://codeload.github.com/", reference[1L],
                               "/tar.gz/", reference[2L]),
                        paste0(name, "_", entry$Version, ".tar.gz"), entry$Archive_SHA256)
    return(paste0(name, "=local::", unpack(archive, name, name)))
  }
  if (startsWith(entry$Source, "Bioconductor")) {
    filename <- paste0(name, "_", entry$Version, ".tar.gz")
    archive <- download(paste0(
      c("https://bioconductor.posit.co", "https://bioconductor.org"),
      "/packages/3.22/bioc/src/contrib/", filename), filename, entry$Archive_SHA256)
  } else {
    archive <- cran_tarball(name, entry$Version, entry$Archive_SHA256)
  }
  paste0(name, "=local::", archive)
}, character(1L))
if (length(refs)) {
  pak::pkg_install(refs, lib = lib, upgrade = FALSE, ask = FALSE, dependencies = NA)
}
writeLines(pin_lines, file.path(find.package("lisi", lib.loc = lib), "BRAINOMICS_PINNED_SOURCE"))
installed <- utils::installed.packages(lib.loc = .libPaths())
installed <- installed[!duplicated(installed[, "Package"]), c("Package", "Version", "LibPath"), drop = FALSE]
utils::write.table(installed, file.path(dirname(lib), "installed-packages.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
unlink(work, recursive = TRUE)
cat("R packages restored.\n")
