#!/usr/bin/env Rscript
# Build the per-cell annotation table used by the export: Cells, Cluster, Dataset,
# CellType and Annotation_Version, in the package cell order.
#
# Usage: export_cell_annotations.R FROZEN_DIR ANALYSIS_DIR WORK_DIR
suppressPackageStartupMessages({library(data.table)})
setDTthreads(4L)
a <- commandArgs(TRUE)
stopifnot(length(a) == 3L)
frozen <- normalizePath(a[1])
analysis <- normalizePath(a[2])
work <- a[3]
dir.create(work, recursive = TRUE, showWarnings = FALSE)

clusters <- as.data.table(readRDS(file.path(frozen, 'cluster_assignments.rds')))
stopifnot(nrow(clusters) == 2602031L, !anyDuplicated(clusters$Cell))
find_one <- function(root, filename) {
  candidates <- list.files(root, recursive = TRUE, full.names = TRUE)
  matches <- candidates[basename(candidates) == filename]
  if (length(matches) != 1L) {
    stop('Expected one ', filename, ' under ', root, '; found ', length(matches))
  }
  matches[[1L]]
}
metadata_file <- Sys.getenv('BRAINOMICS_METADATA_FILE', unset = '')
if (!nzchar(metadata_file)) {
  metadata_file <- find_one(analysis, 'metadata_working.rds')
}
metadata <- as.data.table(readRDS(metadata_file))[, .(Cells, Dataset)]
stopifnot(nrow(metadata) == nrow(clusters), !anyDuplicated(metadata$Cells))
metadata <- metadata[match(clusters$Cell, Cells)]
stopifnot(identical(metadata$Cells, clusters$Cell))

annotation_file <- Sys.getenv('BRAINOMICS_CLUSTER_ANNOTATION', unset = '')
if (!nzchar(annotation_file)) {
  annotation_file <- find_one(analysis, 'final_cluster_annotation.tsv')
}
annotation <- fread(annotation_file, select = c('Cluster', 'CellType', 'Annotation_Version'))
stopifnot(!anyDuplicated(annotation$Cluster))

index <- match(as.character(clusters$Cluster), as.character(annotation$Cluster))
stopifnot(!anyNA(index))
out <- data.table(
  Cells = clusters$Cell,
  Cluster = as.character(clusters$Cluster),
  Dataset = metadata$Dataset,
  CellType = annotation$CellType[index],
  Annotation_Version = annotation$Annotation_Version[index]
)
stopifnot(!anyNA(out), uniqueN(out$CellType) == 12L)
tmp <- file.path(work, 'final_cell_annotations.partial.tsv.gz')
fwrite(out, tmp, sep = '\t', quote = FALSE, na = '')
stopifnot(file.rename(tmp, file.path(work, 'final_cell_annotations.tsv.gz')))
message('Exported ', nrow(out), ' cell annotations')
