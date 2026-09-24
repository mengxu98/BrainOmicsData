#!/usr/bin/env Rscript
# Rebuild the current formal annotation without reading historical cohorts.
suppressPackageStartupMessages(library(data.table))

source_file <- file.path(
  "results", "analysis_run", "07_downstream", "revision_20260918",
  "18_final_annotation_20260918", "final_cell_annotations.tsv.gz"
)
output_dir <- file.path("results", "annotation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation <- fread(source_file, select = c("Cells", "Cluster", "CellType"))
mapping <- fread(file.path(output_dir, "cluster_annotation.tsv"))
stopifnot(
  nrow(annotation) == 2602031L,
  !anyNA(annotation), !anyDuplicated(annotation$Cells),
  uniqueN(annotation$Cluster) == 75L,
  uniqueN(annotation$CellType) == 12L,
  nrow(mapping) == 75L, !anyDuplicated(mapping$Cluster)
)
index <- match(annotation$Cluster, mapping$Cluster)
stopifnot(!anyNA(index),
          identical(annotation$CellType, mapping$CellType[index]),
          identical(as.integer(table(factor(annotation$Cluster,
                                            levels = mapping$Cluster))),
                    as.integer(mapping$Cells)))

assignments_file <- file.path(output_dir, "celltype_assignments.rds")
assignments <- as.data.frame(annotation)
if (!file.exists(assignments_file) ||
    !identical(readRDS(assignments_file), assignments)) {
  temporary <- paste0(assignments_file, ".tmp")
  saveRDS(assignments, temporary, compress = FALSE)
  stopifnot(file.rename(temporary, assignments_file))
}

type_order <- unique(mapping$CellType)
counts <- data.table(
  CellType = type_order,
  Cells = as.integer(table(factor(annotation$CellType, levels = type_order)))
)
stopifnot(nrow(counts) == 12L, sum(counts$Cells) == 2602031L)
fwrite(counts, file.path(output_dir, "celltype_counts.tsv"), sep = "\t")
message("Current formal annotation verified: 2,602,031 cells, 75 clusters, 12 types")
