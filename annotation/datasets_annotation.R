#!/usr/bin/env Rscript
# Prepare summaries for the adopted cell-type annotation.
suppressPackageStartupMessages(library(data.table))

output_dir <- file.path("results", "annotation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

assignments_file <- file.path(output_dir, "celltype_assignments.rds")
if (!file.exists(assignments_file)) {
  stop("Missing adopted cell-type assignments: ", assignments_file)
}
annotation <- as.data.table(readRDS(assignments_file))
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

type_order <- unique(mapping$CellType)
counts <- data.table(
  CellType = type_order,
  Cells = as.integer(table(factor(annotation$CellType, levels = type_order)))
)
stopifnot(nrow(counts) == 12L, sum(counts$Cells) == 2602031L)
fwrite(counts, file.path(output_dir, "celltype_counts.tsv"), sep = "\t")
message("Annotation summaries written: 2,602,031 cells, 75 clusters, 12 types")
