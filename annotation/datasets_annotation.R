#!/usr/bin/env Rscript
# Prepare summaries for the adopted cell-type annotation.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(TRUE)
stopifnot(length(args) == 0L || identical(args, "--build"))
output_dir <- Sys.getenv("BRAINOMICS_ANNOTATION_DIR", unset = file.path("results", "annotation"))
mapping_file <- Sys.getenv("BRAINOMICS_ANNOTATION_TABLE", unset = file.path(output_dir, "cluster_annotation.tsv"))
mapping <- fread(mapping_file)
stopifnot(all(c("Cluster", "Cells", "CellType") %in% names(mapping)),
  nrow(mapping) == 75L, !anyDuplicated(mapping$Cluster),
  sum(mapping$Cells) == 2602031L, uniqueN(mapping$CellType) == 12L,
  !anyNA(mapping[, .(Cluster, Cells, CellType)]))

assignments_file <- file.path(output_dir, "celltype_assignments.rds")
if (identical(args, "--build")) {
  if (file.exists(assignments_file)) stop("Output already exists: ", assignments_file)
  cluster_file <- Sys.getenv("BRAINOMICS_CLUSTER_ASSIGNMENTS", unset = "")
  if (!nzchar(cluster_file)) stop("Set BRAINOMICS_CLUSTER_ASSIGNMENTS to the retained adopted partition")
  clusters <- as.data.table(readRDS(cluster_file))
  stopifnot(all(c("Cell", "Cluster") %in% names(clusters)))
  index <- match(as.character(clusters$Cluster), as.character(mapping$Cluster))
  stopifnot(!anyNA(index))
  annotation <- data.table(Cells = as.character(clusters$Cell),
    Cluster = as.character(clusters$Cluster), CellType = mapping$CellType[index])
} else {
  if (!file.exists(assignments_file)) stop("Missing adopted cell-type assignments: ", assignments_file)
  annotation <- as.data.table(readRDS(assignments_file))
}
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
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (identical(args, "--build")) {
  temporary <- paste0(assignments_file, ".tmp.", Sys.getpid())
  saveRDS(as.data.frame(annotation), temporary)
  stopifnot(file.rename(temporary, assignments_file))
}
fwrite(counts, file.path(output_dir, "celltype_counts.tsv"), sep = "\t")
message("Annotation summaries written: 2,602,031 cells, 75 clusters, 12 types")
