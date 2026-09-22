#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

repo_dir <- normalizePath(value_after("--repo-dir", "."), mustWork = TRUE)
integration_dir <- normalizePath(
  value_after("--integration-dir", "../../data/BrainOmicsData/integration"),
  mustWork = TRUE
)
source(file.path(repo_dir, "functions", "metadata_schema.R"))
source(file.path(repo_dir, "functions", "sample_schema.R"))

targets_arg <- value_after(
  "--targets",
  "metadata_filtered.rds,objects_celltype_plot.rds"
)
targets <- trimws(strsplit(targets_arg, ",", fixed = TRUE)[[1]])
backup_dir_arg <- value_after("--backup-dir", "")
backup_dir <- if (nzchar(backup_dir_arg)) normalizePath(backup_dir_arg, mustWork = FALSE) else ""
no_backup <- "--no-backup" %in% args
timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")

save_with_backup <- function(object, path) {
  if (!file.exists(path)) stop("missing RDS file: ", path)
  backup <- NA_character_
  if (!no_backup) {
    backup_base <- paste0(basename(path), ".bak_", timestamp)
    backup <- if (nzchar(backup_dir)) file.path(backup_dir, backup_base) else file.path(dirname(path), backup_base)
    dir.create(dirname(backup), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(path, backup)) {
      stop("failed to move original file to backup: ", path)
    }
  }
  ok <- FALSE
  tryCatch(
    {
      saveRDS(object, path)
      ok <<- TRUE
    },
    error = function(e) {
      if (file.exists(path)) unlink(path)
      if (!no_backup && file.exists(backup)) file.rename(backup, path)
      stop(e)
    }
  )
  if (ok) {
    message("wrote: ", path)
    if (!no_backup) {
      message("backup: ", backup)
    } else {
      message("backup: disabled")
    }
  }
}

update_rds <- function(path) {
  message("loading: ", path)
  object <- readRDS(path)
  if (is.data.frame(object)) {
    object <- add_metadata_schema(object)
    object <- add_sample_schema(object)
  } else if (inherits(object, "Seurat")) {
    object@meta.data <- add_metadata_schema(object@meta.data)
    object@meta.data <- add_sample_schema(object@meta.data)
  } else {
    stop("unsupported RDS class for ", path, ": ", paste(class(object), collapse = ", "))
  }
  save_with_backup(object, path)
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("The Seurat package is required to update Seurat RDS metadata.")
}

files <- file.path(integration_dir, targets)
for (file in files) update_rds(file)
