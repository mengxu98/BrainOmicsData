#!/usr/bin/env Rscript
# Seed and resolution sensitivity analysis of the adopted clustering.
#
# The analysis re-runs the documented clustering procedure on the 2,602,031-cell
# cohort and reports the cluster count, the ARI against the adopted 75-cluster
# partition and the per-cell-type majority agreement.
#
# Design (see the clustering contract for the same parameters):
#   cohort       2,602,031 cells
#   graph        symmetric adaptive Gaussian weighted kNN, 20 neighbours, 50 RPCA dimensions
#   algorithm    Louvain (Seurat algorithm 1), n.start = 10, n.iter = 10, group.singletons = TRUE
#   runs         seed 20260730 / resolution 2   (baseline partition)
#                seed 42       / resolution 2
#                seed 2026     / resolution 2
#                seed 20260730 / resolution 1.5
#                seed 20260730 / resolution 2.5
#   metrics      cluster count; ARI against the 75-cluster partition
#                per-cell-type majority mapped agreement against the 12-type annotation
#
# Published values to check against:
#   75 clusters ARI 1        (baseline replay)
#   77 clusters ARI 0.754474 (seed 42)
#   72 clusters ARI 0.803025 (seed 2026)
#   64 clusters ARI 0.712368 (resolution 1.5)
#   86 clusters ARI 0.685946 (resolution 2.5)
#
# Inputs
#   --frozen DIR   directory holding cluster_assignments.rds and cell_metadata.rds
#                  (default: the frozen results directory of the current release)
#   --partitions DIR   optional directory holding precomputed partitions as
#                  <seed>_<resolution>.rds, each a data.table with Cells and Cluster.
#                  One of --partitions or --graph is required for the four
#                  non-baseline runs.
#   --graph FILE   optional Seurat Graph object, identical to the production graph.
#                  When supplied, the four non-baseline partitions are recomputed.
#   --out DIR      output directory (default: <frozen>/seed_resolution_sensitivity)
#   --annotation FILE   adopted cell-type annotation table with Cells, Cluster, CellType
#                  (default: derived from the frozen metadata plus
#                  final_cluster_annotation.tsv when present)
#
# Usage
#   Rscript reconstruct_seed_resolution_sensitivity.R --frozen <dir> [--graph <graph.rds>] --out <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
})

setDTthreads(2L)

args <- commandArgs(trailingOnly = TRUE)
value_of <- function(flag, default = NA_character_) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  stopifnot(hit < length(args))
  args[[hit + 1L]]
}

frozen <- normalizePath(value_of("--frozen"), mustWork = TRUE)
graph_file <- value_of("--graph")
partitions_dir <- value_of("--partitions")
annotation_file <- value_of("--annotation")
out <- value_of("--out", file.path(frozen, "seed_resolution_sensitivity"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

runs <- data.table(
  Task = 1:5,
  Seed = c(20260730L, 42L, 2026L, 20260730L, 20260730L),
  Resolution = c(2, 2, 2, 1.5, 2.5),
  Role = c("baseline_replay", "seed", "seed", "resolution", "resolution"),
  Published_Clusters = c(75L, 77L, 72L, 64L, 86L),
  Published_ARI = c(1, 0.754473754912277, 0.803025492214823, 0.712367845263224, 0.685946143287341)
)

adjusted_rand_index <- function(left, right) {
  stopifnot(length(left) == length(right), length(left) > 1L, !anyNA(left), !anyNA(right))
  choose2 <- function(x) x * (x - 1) / 2
  tab <- table(left, right)
  observed <- sum(choose2(tab))
  rows <- sum(choose2(rowSums(tab)))
  cols <- sum(choose2(colSums(tab)))
  expected <- rows * cols / choose2(sum(tab))
  maximum <- (rows + cols) / 2
  if (maximum == expected) 1 else (observed - expected) / (maximum - expected)
}

majority_agreement <- function(reference, labels) {
  stopifnot(length(reference) == length(labels), length(reference) > 1L)
  counts <- data.table(Type = as.character(reference), Cluster = as.character(labels))[
    , .(Overlap = .N), by = .(Type, Cluster)]
  mapping <- counts[order(-Overlap, Cluster)][, .SD[1L], by = Cluster]
  predicted <- mapping$Type[match(labels, mapping$Cluster)]
  dt <- data.table(CellType = as.character(reference), Predicted = predicted)
  per_type <- dt[, .(Cells = .N, Agreement = mean(Predicted == CellType)), by = CellType]
  weighted <- dt[, mean(Predicted == CellType)]
  list(weighted = weighted, per_type = per_type)
}

clusters <- as.data.table(readRDS(file.path(frozen, "cluster_assignments.rds")))
stopifnot(all(c("Cell", "Cluster") %in% names(clusters)), nrow(clusters) == 2602031L)
clusters[, Cluster := as.character(Cluster)]

read_cluster_annotation <- function(path) {
  dt <- fread(path)
  stopifnot(all(c("Cluster", "CellType") %in% names(dt)))
  unique(dt[, .(Cluster = as.character(Cluster), CellType = as.character(CellType))])
}

if (!is.na(annotation_file)) {
  cluster_annotation <- read_cluster_annotation(annotation_file)
} else {
  candidates <- c(
    file.path(frozen, "final_cluster_annotation.tsv"),
    file.path(frozen, "cluster_annotation.tsv"),
    file.path(frozen, "annotation", "final_cluster_annotation.tsv"),
    file.path(dirname(frozen), "analysis_from_scratch_20260917", "07_downstream",
      "revision_20260918", "18_final_annotation_20260918", "final_cluster_annotation.tsv"),
    file.path(dirname(frozen), "analysis_from_scratch_20260917", "script_output",
      "cluster_annotation.tsv")
  )
  # Fall back to any annotation table with Cluster and CellType under the sibling
  # analysis directory, so the script works without hardcoded deep paths.
  analysis_root <- file.path(dirname(frozen), "analysis_from_scratch_20260917")
  if (dir.exists(analysis_root)) {
    discovered <- list.files(analysis_root, pattern = "cluster_annotation[.]tsv$",
      recursive = TRUE, full.names = TRUE)
    candidates <- unique(c(candidates, discovered))
  }
  annotation_path <- candidates[file.exists(candidates)][1L]
  if (is.na(annotation_path)) {
    stop("No cluster annotation table with Cluster and CellType was found. ",
      "Pass --annotation <file>.")
  }
  cluster_annotation <- read_cluster_annotation(annotation_path)
}

stopifnot(uniqueN(cluster_annotation$Cluster) == 75L, uniqueN(cluster_annotation$CellType) == 12L)
annotation <- cluster_annotation[match(clusters$Cluster, Cluster)]
stopifnot(!anyNA(annotation$CellType))
annotation <- data.table(CellType = annotation$CellType)
stopifnot(nrow(annotation) == nrow(clusters))

graph <- NULL
if (!is.na(graph_file)) {
  graph <- readRDS(graph_file)
  stopifnot(inherits(graph, "Graph"), identical(rownames(graph), clusters$Cell))
}

read_partition <- function(seed, resolution) {
  if (seed == 20260730L && resolution == 2) {
    return(list(Cluster = clusters$Cluster, origin = "published frozen partition"))
  }
  if (!is.na(partitions_dir)) {
    path <- file.path(partitions_dir, sprintf("%d_%s.rds", seed, resolution))
    if (file.exists(path)) {
      part <- as.data.table(readRDS(path))
      stopifnot(all(c("Cells", "Cluster") %in% names(part)))
      part <- part[match(clusters$Cell, Cells)]
      stopifnot(identical(as.character(part$Cells), as.character(clusters$Cell)))
      return(list(Cluster = as.character(part$Cluster), origin = path))
    }
  }
  if (!is.null(graph)) {
    set.seed(seed)
    partition <- FindClusters(graph, algorithm = 1L, resolution = resolution,
      n.start = 10L, n.iter = 10L, random.seed = seed, group.singletons = TRUE,
      verbose = FALSE)
    return(list(Cluster = as.character(partition[[1L]]), origin = "recomputed from supplied graph"))
  }
  stop(sprintf(paste0("No partition for seed %d / resolution %s. Supply --partitions ",
    "with the stored partitions or --graph to recompute the partition."),
    seed, format(resolution)))
}

summary_rows <- list()
type_rows <- list()
for (i in seq_len(nrow(runs))) {
  run <- runs[i]
  part <- read_partition(run$Seed, run$Resolution)
  labels <- part$Cluster
  stopifnot(length(labels) == nrow(clusters), !anyNA(labels))
  ari <- adjusted_rand_index(clusters$Cluster, labels)
  agreement <- majority_agreement(annotation$CellType, labels)
  summary_rows[[i]] <- data.table(
    Task = run$Task, Seed = run$Seed, Resolution = run$Resolution, Role = run$Role,
    Cells = nrow(clusters), Clusters = uniqueN(labels),
    ARI_To_Published_Partition = ari,
    CellType_Majority_Agreement = agreement$weighted,
    Published_Clusters = run$Published_Clusters, Published_ARI = run$Published_ARI,
    Origin = part$origin
  )
  type_rows[[i]] <- data.table(Task = run$Task, Seed = run$Seed, Resolution = run$Resolution,
    agreement$per_type)
}

summary <- rbindlist(summary_rows)
summary[, Clusters_Match_Published := Clusters == Published_Clusters]
summary[, ARI_Delta_From_Published := ARI_To_Published_Partition - Published_ARI]
by_type <- rbindlist(type_rows)

fwrite(summary, file.path(out, "stability_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(by_type, file.path(out, "stability_by_celltype.tsv"), sep = "\t", quote = FALSE)
fwrite(runs, file.path(out, "stability_design.tsv"), sep = "\t", quote = FALSE)
writeLines(c(
  paste0("Status=complete_or_partial; generated=", format(Sys.time())),
  "Script_Status=reconstructed_from_documented_parameters; the original script was not archived",
  paste0("Cells=", nrow(clusters)),
  paste0("Runs_Written=", nrow(summary)),
  "Reference=Supplementary_Data.xlsx sheets Clustering_Stability and Stability_By_Cell_Type"
), file.path(out, "run_info.txt"))

print(summary)
