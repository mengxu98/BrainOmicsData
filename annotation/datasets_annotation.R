#!/usr/bin/env Rscript
# Apply cluster labels without changing expression or embeddings.
suppressPackageStartupMessages(library(SeuratObject))
source("functions/utils.R")
metadata <- readRDS("../../data/BrainOmicsData/integration_25/evaluation/plot_metadata_slim.rds")
mapping <- read_tsv("results/annotation/cluster_annotation.tsv")
index <- match(metadata$Cluster, mapping$Cluster)
stopifnot(
  nrow(mapping) == 75L, !anyDuplicated(mapping$Cluster), !anyNA(index),
  nrow(metadata) == 2602031L, !anyDuplicated(metadata$Cells)
)
assignments <- data.frame(
  Cells = as.character(metadata$Cells),
  Cluster = as.character(metadata$Cluster), CellType = mapping$CellType[index]
)
path <- "../../data/BrainOmicsData/integration_25/annotation/celltype_assignments.rds"
if (file.exists(path)) {
  current <- as.data.frame(readRDS(path))
  stopifnot(identical(current$Cells, assignments$Cells))
} else {
  current <- NULL
}
if (!identical(current, assignments)) {
  saveRDS(assignments, paste0(path, ".tmp"), compress = FALSE)
  stopifnot(file.rename(paste0(path, ".tmp"), path))
}
type_order <- unique(mapping$CellType)
counts <- data.frame(
  CellType = type_order,
  Cells = as.integer(table(factor(assignments$CellType, levels = type_order)))
)
stopifnot(
  nrow(counts) == 12L, sum(counts$Cells) == 2602031L,
  counts$Cells[counts$CellType == "Excitatory neurons"] == 971364L
)
for (output in c(
  "results/annotation/celltype_counts.tsv",
  "../../data/BrainOmicsData/integration_25/annotation/celltype_counts.tsv"
)) {
  write_tsv(counts, output)
}
rm(metadata, current)
gc()
for (name in c(
  "metadata_filtered.rds", "evaluation/plot_metadata_slim.rds",
  "objects_celltype_plot.rds"
)) {
  path <- file.path("../../data/BrainOmicsData/integration_25", name)
  object <- readRDS(path)
  seurat <- inherits(object, "Seurat")
  meta <- if (seurat) object[[]] else as.data.frame(object)
  cells <- if (seurat) colnames(object) else as.character(meta$Cells)
  index <- match(cells, assignments$Cells)
  stopifnot(length(cells) == nrow(assignments), !anyNA(index), !anyDuplicated(cells))
  fields <- c("Cluster", "CellType")
  changed <- any(vapply(fields, function(field) {
    !identical(as.character(meta[[field]]), assignments[[field]][index])
  }, logical(1)))
  obsolete <- intersect(c("Detailed_CellType", "Main_CellType", "CellType_Broad"), names(meta))
  changed <- changed || length(obsolete) > 0L
  if (changed) {
    for (field in fields) meta[[field]] <- assignments[[field]][index]
    for (field in obsolete) meta[[field]] <- NULL
    if (seurat) {
      object@meta.data <- meta
    } else {
      if (inherits(object, "data.table")) meta <- data.table::as.data.table(meta)
      object <- meta
    }
    validate_celltype_metadata(meta, cells)
    saveRDS(object, paste0(path, ".tmp"), compress = FALSE)
    stopifnot(file.rename(paste0(path, ".tmp"), path))
  }
  message(name, if (changed) ": labels updated" else ": labels already match")
  rm(object, meta, cells, index)
  gc()
}
