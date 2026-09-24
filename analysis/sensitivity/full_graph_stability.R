#!/usr/bin/env Rscript
# Louvain seed and resolution sensitivity on the fixed full-cohort RPCA graph.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})
source("functions/data_paths.R")

args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
input_dir <- normalizePath(brainomics_data_path(args[1]), mustWork = TRUE)
output_dir <- normalizePath(args[2])
task <- as.integer(args[3])

settings <- data.table(
  Task = 1:5,
  Seed = c(20260730L, 42L, 2026L, 20260730L, 20260730L),
  Resolution = c(2, 2, 2, 1.5, 2.5),
  Role = c("reference", "seed", "seed", "resolution", "resolution")
)
stopifnot(task %in% settings$Task)
setting <- settings[Task == task]
run_dir <- file.path(output_dir, paste0("run_", task))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

if (task > 1L) {
  reference_summary <- fread(file.path(output_dir, "run_1", "summary.tsv"))
  stopifnot(nrow(reference_summary) == 1L, reference_summary$Task == 1L,
            abs(reference_summary$ARI - 1) < 1e-12)
}

graph_file <- file.path(input_dir, "rpca_symmetric_weighted_knn_graph.rds")
partition_file <- file.path(input_dir, "cluster_assignments.rds")
reference <- as.data.table(readRDS(partition_file))
graph <- readRDS(graph_file)
stopifnot(
  nrow(reference) == 2602031L,
  uniqueN(reference$Cluster) == 75L,
  !anyDuplicated(reference$Cell),
  is(graph, "Graph"),
  identical(rownames(graph), reference$Cell),
  identical(colnames(graph), reference$Cell),
  isSymmetric(graph), all(is.finite(graph@x)), all(graph@x > 0)
)

parameters <- data.table(
  Task = task, Seed = setting$Seed, Resolution = setting$Resolution,
  Role = setting$Role, Cells = nrow(reference),
  Graph_MD5 = unname(tools::md5sum(graph_file)),
  Reference_Assignment_MD5 = unname(tools::md5sum(partition_file)),
  Algorithm = 1L, Starts = 10L, Iterations = 10L,
  Group_Singletons = TRUE, Subsampling = FALSE
)
fwrite(parameters, file.path(run_dir, "parameters.tsv"), sep = "\t")

future::plan("sequential")
set.seed(setting$Seed)
started <- proc.time()[["elapsed"]]
partition <- FindClusters(
  graph, algorithm = 1L, resolution = setting$Resolution,
  n.start = 10L, n.iter = 10L, random.seed = setting$Seed,
  group.singletons = TRUE, verbose = TRUE
)
seconds <- proc.time()[["elapsed"]] - started
labels <- paste0("C", as.character(partition[[1L]]))
stopifnot(length(labels) == nrow(reference), !anyNA(labels))
rm(graph, partition)
gc()

overlap <- table(Baseline = reference$Cluster, Perturbed = labels)
choose2 <- function(value) value * (value - 1) / 2
cells <- sum(overlap)
row_pairs <- sum(choose2(rowSums(overlap)))
column_pairs <- sum(choose2(colSums(overlap)))
observed_pairs <- sum(choose2(overlap))
expected_pairs <- row_pairs * column_pairs / choose2(cells)
ari <- (observed_pairs - expected_pairs) /
  ((row_pairs + column_pairs) / 2 - expected_pairs)

summary <- data.table(
  Task = task, Seed = setting$Seed, Resolution = setting$Resolution,
  Role = setting$Role, Cells = cells, Clusters = uniqueN(labels), ARI = ari,
  Seconds = seconds
)
if (task == 1L) {
  stopifnot(abs(ari - 1) < 1e-12)
}
fwrite(summary, file.path(run_dir, "summary.tsv"), sep = "\t")
overlap_table <- as.data.table(overlap)
overlap_table <- overlap_table[N > 0L]
setnames(overlap_table, "N", "Freq")
fwrite(overlap_table, file.path(run_dir, "overlap.tsv"), sep = "\t")

cluster_metrics <- rbindlist(lapply(rownames(overlap), function(cluster) {
  values <- overlap[cluster, ]
  size <- sum(values)
  data.table(
    Baseline = cluster, Cells = size,
    Largest_Fragment_Fraction = max(values) / size,
    Best_Jaccard = max(values / (size + colSums(overlap) - values)),
    Fragments = sum(values > 0)
  )
}))
fwrite(cluster_metrics, file.path(run_dir, "cluster_stability.tsv"), sep = "\t")
saveRDS(data.table(Cell = reference$Cell, Cluster = labels),
        file.path(run_dir, "assignments.rds"), compress = "gzip")
write_json(
  list(state = "COMPLETE", cells = cells, task = task, ARI = round(ari, 4L)),
  file.path(run_dir, "COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE
)
message("Clustering sensitivity result written for task ", task)
