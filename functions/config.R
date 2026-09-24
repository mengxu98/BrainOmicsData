# Shared configuration for the manuscript figure scripts.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
setDTthreads(2L)
source("functions/utils.R")

root <- normalizePath(
  Sys.getenv("BRAINOMICS_ANALYSIS_DIR", unset = file.path("results", "analysis_run")),
  mustWork = TRUE
)

resolve_directory <- function(variable, candidates, marker, label) {
  configured <- Sys.getenv(variable, unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))
  candidates <- unique(candidates[dir.exists(candidates)])
  matches <- candidates[file.exists(file.path(candidates, marker))]
  if (length(matches) != 1L) {
    stop("Set ", variable, " to the directory containing ", marker,
         " (found ", length(matches), " ", label, " directories)")
  }
  normalizePath(matches[[1L]], mustWork = TRUE)
}

resolve_input_file <- function(variable, search_root, filename) {
  configured <- Sys.getenv(variable, unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))
  candidates <- list.files(search_root, recursive = TRUE, full.names = TRUE)
  matches <- candidates[basename(candidates) == filename]
  if (length(matches) != 1L) {
    stop("Set ", variable, " to ", filename,
         " (found ", length(matches), " matching files)")
  }
  normalizePath(matches[[1L]], mustWork = TRUE)
}

downstream_root <- file.path(root, "07_downstream")
run <- resolve_directory(
  "BRAINOMICS_ANALYSIS_RUN_DIR",
  c(root, if (dir.exists(downstream_root)) list.dirs(downstream_root, recursive = FALSE, full.names = TRUE)),
  file.path("01_metadata", "metadata_working.rds"),
  "analysis"
)
doc <- resolve_directory(
  "BRAINOMICS_FIGURE_DATA_DIR",
  c(run, list.dirs(run, recursive = FALSE, full.names = TRUE)),
  file.path("tables", "age_signal", "age_metrics_summary.tsv"),
  "figure-data"
)

annotation_file <- Sys.getenv(
  "BRAINOMICS_ANNOTATION_TABLE",
  unset = file.path("results", "annotation", "cluster_annotation.tsv")
)
a <- fread(annotation_file)
stopifnot(all(c("Cluster", "Cells", "CellType") %in% names(a)))
a[, Working_CellType := CellType]
a[, Colour := unname(brainomics_celltype_colors[Working_CellType])]
type_levels <- names(brainomics_celltype_colors)
stopifnot(
  nrow(a) == 75L,
  sum(a$Cells) == 2602031L,
  !anyNA(brainomics_celltype_colors)
)
