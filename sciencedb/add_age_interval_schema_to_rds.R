#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

integration_dir <- normalizePath(
  value_after("--integration-dir", "../../data/BrainOmicsData/integration"),
  mustWork = TRUE
)
targets_arg <- value_after(
  "--targets",
  "metadata_filtered.rds,objects_celltypes.rds,objects_celltype_plot.rds"
)
targets <- trimws(strsplit(targets_arg, ",", fixed = TRUE)[[1]])
timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")

age_interval <- c(
  "S1" = "Embryonic",
  "S2" = "Early fetal",
  "S3" = "Early fetal",
  "S4" = "Early mid-fetal",
  "S5" = "Early mid-fetal",
  "S6" = "Late mid-fetal",
  "S7" = "Late fetal",
  "S8" = "Neonatal and early infancy",
  "S9" = "Late infancy",
  "S10" = "Early childhood",
  "S11" = "Middle and late childhood",
  "S12" = "Adolescence",
  "S13" = "Young adulthood",
  "S14" = "Middle adulthood",
  "S15" = "Late adulthood"
)
age_range <- c(
  "S1" = "4-8 PCW",
  "S2" = "8-10 PCW",
  "S3" = "10-13 PCW",
  "S4" = "13-16 PCW",
  "S5" = "16-19 PCW",
  "S6" = "19-24 PCW",
  "S7" = "24-38 PCW",
  "S8" = "0-0.5 years",
  "S9" = "0.5-1 years",
  "S10" = "1-6 years",
  "S11" = "6-12 years",
  "S12" = "12-20 years",
  "S13" = "20-40 years",
  "S14" = "40-60 years",
  "S15" = "60+ years"
)

add_age_interval_schema <- function(meta) {
  if (!"Stage" %in% names(meta)) {
    stop("metadata is missing required legacy column: Stage")
  }
  normalize_dashes <- function(x) {
    if (!is.character(x)) return(x)
    x <- gsub(intToUtf8(0x2013), "-", x, fixed = TRUE)
    x <- gsub(intToUtf8(0x2014), "-", x, fixed = TRUE)
    x <- gsub(intToUtf8(0x2212), "-", x, fixed = TRUE)
    x
  }
  meta[] <- lapply(meta, normalize_dashes)
  meta$AgeIntervalID <- as.character(meta$Stage)
  meta$AgeInterval <- unname(age_interval[meta$AgeIntervalID])
  meta$AgeRange <- unname(age_range[meta$AgeIntervalID])
  missing_interval <- is.na(meta$AgeInterval) | is.na(meta$AgeRange)
  if (any(missing_interval)) {
    bad <- unique(meta$AgeIntervalID[missing_interval])
    stop("unmapped Stage values: ", paste(bad, collapse = ", "))
  }
  meta
}

save_with_backup <- function(object, path) {
  backup <- paste0(path, ".bak_", timestamp)
  if (!file.exists(path)) stop("missing RDS file: ", path)
  if (!file.rename(path, backup)) {
    stop("failed to move original file to backup: ", path)
  }
  ok <- FALSE
  tryCatch(
    {
      saveRDS(object, path)
      ok <<- TRUE
    },
    error = function(e) {
      if (file.exists(path)) unlink(path)
      file.rename(backup, path)
      stop(e)
    }
  )
  if (ok) {
    message("wrote: ", path)
    message("backup: ", backup)
  }
}

update_rds <- function(path) {
  message("loading: ", path)
  object <- readRDS(path)
  if (is.data.frame(object)) {
    object <- add_age_interval_schema(object)
  } else if (inherits(object, "Seurat")) {
    object@meta.data <- add_age_interval_schema(object@meta.data)
  } else {
    stop("unsupported RDS class for ", path, ": ", paste(class(object), collapse = ", "))
  }
  save_with_backup(object, path)
}

files <- file.path(
  integration_dir,
  targets
)
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("The Seurat package is required to update Seurat RDS metadata.")
}
for (file in files) update_rds(file)
