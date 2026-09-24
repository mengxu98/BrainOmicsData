#!/usr/bin/env Rscript
# Summarize the five clustering-sensitivity partitions.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
source("functions/data_paths.R")

args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
input_dir <- normalizePath(brainomics_data_path(args[1]), mustWork = TRUE)
output_dir <- normalizePath(args[2])

reference <- as.data.table(readRDS(file.path(input_dir, "cluster_assignments.rds")))
annotation <- fread(file.path(output_dir, "final_cluster_annotation.tsv"))
stopifnot(nrow(reference) == 2602031L, nrow(annotation) == 75L,
          uniqueN(annotation$CellType) == 12L)
dataset <- sub("::.*$", "", reference$Cell)
stopifnot(uniqueN(dataset) == 22L)
cell_type <- annotation$CellType[match(reference$Cluster, annotation$Cluster)]
stopifnot(!anyNA(cell_type))

adjusted_rand_index <- function(left, right) {
  tab <- table(left, right)
  choose2 <- function(value) value * (value - 1) / 2
  cells <- sum(tab)
  row_pairs <- sum(choose2(rowSums(tab)))
  column_pairs <- sum(choose2(colSums(tab)))
  expected <- row_pairs * column_pairs / choose2(cells)
  denominator <- (row_pairs + column_pairs) / 2 - expected
  if (denominator == 0) return(if (identical(left, right)) 1 else NA_real_)
  (sum(choose2(tab)) - expected) / denominator
}

summaries <- list()
by_dataset <- list()
by_celltype <- list()
by_cluster <- list()
for (task in 1:5) {
  run_dir <- file.path(output_dir, paste0("run_", task))
  completion <- fromJSON(file.path(run_dir, "COMPLETE.json"))
  stopifnot(
    identical(completion$state, "COMPLETE"),
    completion$task == task,
    completion$cells == nrow(reference)
  )
  assignments <- as.data.table(readRDS(file.path(run_dir, "assignments.rds")))
  stopifnot(identical(assignments$Cell, reference$Cell))
  summaries[[task]] <- fread(file.path(run_dir, "summary.tsv"))
  stopifnot(abs(completion$ARI - round(summaries[[task]]$ARI, 4L)) < 1e-12)

  counts <- table(cell_type, assignments$Cluster)
  majority <- rownames(counts)[max.col(t(counts), ties.method = "first")]
  names(majority) <- colnames(counts)
  mapped <- unname(majority[assignments$Cluster])

  by_celltype[[task]] <- rbindlist(lapply(sort(unique(cell_type)), function(type) {
    selected <- cell_type == type
    data.table(Task = task, CellType = type, Cells = sum(selected),
               Majority_Mapped_Agreement = mean(mapped[selected] == cell_type[selected]))
  }))
  by_dataset[[task]] <- rbindlist(lapply(sort(unique(dataset)), function(source) {
    selected <- dataset == source
    data.table(
      Task = task, Dataset = source, Cells = sum(selected),
      Partition_ARI = adjusted_rand_index(reference$Cluster[selected],
                                           assignments$Cluster[selected]),
      Majority_Mapped_Class_Agreement = mean(mapped[selected] == cell_type[selected])
    )
  }))
  cluster_table <- fread(file.path(run_dir, "cluster_stability.tsv"))
  cluster_table[, Task := task]
  by_cluster[[task]] <- cluster_table
}

fwrite(rbindlist(summaries), file.path(output_dir, "stability_summary.tsv"), sep = "\t")
fwrite(rbindlist(by_dataset), file.path(output_dir, "stability_by_dataset.tsv"), sep = "\t")
fwrite(rbindlist(by_celltype), file.path(output_dir, "stability_by_celltype.tsv"), sep = "\t")
fwrite(rbindlist(by_cluster), file.path(output_dir, "stability_by_cluster.tsv"), sep = "\t")
write_json(
  list(
    state = "PASS", runs = 5L, cells = nrow(reference),
    datasets = uniqueN(dataset), baseline_replay = TRUE,
    scope = paste(
      "All cells on one fixed production graph.",
      "Majority label matching measures stability conditional on the adopted annotation",
      "and is not independent biological accuracy.",
      "Integration randomness and graph construction were not perturbed."
    )
  ),
  file.path(output_dir, "SUMMARY_COMPLETE.json"),
  pretty = TRUE, auto_unbox = TRUE
)
message("Clustering sensitivity summaries written")
