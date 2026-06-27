#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
script_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()
source_file <- file.path(script_dir, "lightweight_scop_h5ad_functions.R")
if (!file.exists(source_file)) {
  stop("Missing lightweight scop source file: ", source_file)
}

pkg_dir <- file.path(tempdir(), "scop")
unlink(pkg_dir, recursive = TRUE, force = TRUE)
dir.create(file.path(pkg_dir, "R"), recursive = TRUE, showWarnings = FALSE)

writeLines(
  c(
    "Package: scop",
    "Type: Package",
    "Title: Lightweight h5ad Export Helpers from scop",
    "Version: 0.8.9.9000",
    "Authors@R: person('Mengxu', 'Xu', role = c('aut', 'cre'), email = 'mengxu@example.com')",
    "Description: Minimal h5ad export subset for BrainOmicsData HPC conversion.",
    "License: MIT",
    "Encoding: UTF-8",
    "Imports: Matrix, reticulate, SeuratObject",
    "NeedsCompilation: no"
  ),
  file.path(pkg_dir, "DESCRIPTION")
)

writeLines(
  c(
    "export(srt_to_adata)",
    "export(srt_to_h5ad)"
  ),
  file.path(pkg_dir, "NAMESPACE")
)

file.copy(source_file, file.path(pkg_dir, "R", "h5ad.R"), overwrite = TRUE)
install.packages(pkg_dir, repos = NULL, type = "source")
