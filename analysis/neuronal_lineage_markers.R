#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(matrixStats)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/processed_object.R")
source("functions/dataset_metadata.R")
source("functions/integration.R")

lineage_id <- Sys.getenv("BRAINOMICS_LINEAGE_ID", unset = "")
lineage_dir <- brainomics_lineage_path(lineage_id)
audit_file <- file.path(
  brainomics_data_path("integration_25"), "lineage_analysis_20260825_v2",
  "lineage_selection.tsv"
)
if (!dir.exists(lineage_dir) || !nzchar(lineage_id)) {
  stop("BRAINOMICS_LINEAGE_ID is required")
}
object_file <- file.path(lineage_dir, "objects_lineage_integrated.rds")
assignment_file <- file.path(
  lineage_dir, "clustering", "lineage_cluster_assignments.tsv.gz"
)
clustering_manifest_file <- file.path(
  lineage_dir, "clustering", "lineage_clustering_output_manifest.tsv"
)
required_files <- c(
  object_file, assignment_file, clustering_manifest_file, audit_file
)
if (any(!file.exists(required_files))) {
  stop(
    "Lineage marker inputs are missing: ",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
}
output_dir <- file.path(lineage_dir, "markers")
checkpoint_dir <- file.path(output_dir, "dataset_checkpoints")
existing_outputs <- if (dir.exists(output_dir)) {
  setdiff(list.files(output_dir, all.files = FALSE), "dataset_checkpoints")
} else {
  character()
}
if (length(existing_outputs) > 0L) {
  stop("Lineage marker workflow refuses to overwrite existing outputs")
}
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

normalization_scale <- 10000
minimum_cluster_cells_per_dataset <- 20L
top_markers_per_direction <- 50L
assignments <- fread(assignment_file, sep = "\t")
required_assignment <- c(
  "Cells", "Global_Cluster", "Main_CellType", "Lineage_Cluster"
)
if (any(!required_assignment %in% names(assignments)) ||
  anyNA(assignments[, ..required_assignment]) ||
  anyDuplicated(assignments$Cells)) {
  stop("Lineage cluster assignment schema changed")
}
main_celltype_audit <- fread(audit_file, sep = "\t")
required_audit_columns <- c(
  "Cluster", "Main_CellType", "Exclude_From_Subtype_Inference"
)
if (any(!required_audit_columns %in% names(main_celltype_audit)) ||
  anyDuplicated(main_celltype_audit$Cluster) ||
  anyNA(main_celltype_audit[, ..required_audit_columns]) ||
  !is.logical(main_celltype_audit$Exclude_From_Subtype_Inference)) {
  stop("Main-cell-type audit failed the subtype-inference contract")
}
lineage_celltype <- unique(as.character(assignments$Main_CellType))
if (length(lineage_celltype) != 1L || is.na(lineage_celltype)) {
  stop("Lineage assignments must contain exactly one main cell type")
}
excluded_subtype_clusters <- main_celltype_audit[
  Main_CellType == lineage_celltype &
    Exclude_From_Subtype_Inference == TRUE,
  as.character(Cluster)
]
unknown_excluded_clusters <- setdiff(
  excluded_subtype_clusters, unique(as.character(assignments$Global_Cluster))
)
if (length(unknown_excluded_clusters) > 0L) {
  stop(
    "Subtype-inference exclusions are absent from the lineage: ",
    paste(unknown_excluded_clusters, collapse = ",")
  )
}
assignments[, Subtype_Inference_Eligible :=
  !Global_Cluster %chin% excluded_subtype_clusters]
marker_assignments <- assignments[Subtype_Inference_Eligible == TRUE]
if (nrow(marker_assignments) == 0L) {
  stop("No lineage cells are eligible for subtype marker inference")
}
cluster_levels <- sort(unique(as.character(marker_assignments$Lineage_Cluster)))
if (length(cluster_levels) < 2L) {
  stop("Lineage marker workflow requires at least two lineage clusters")
}

object <- load_processed_object(object_file)
if (!identical(colnames(object), assignments$Cells) ||
  ncol(object) != nrow(assignments)) {
  stop("Lineage integrated object and assignments differ")
}
features <- rownames(object)
if (length(features) < 10000L || anyNA(features) || anyDuplicated(features)) {
  stop("Lineage object feature space is incomplete")
}
count_layers <- Layers(object, assay = "RNA", search = "^counts")
if (length(count_layers) < 2L) {
  stop("Lineage object does not retain source-dataset count layers")
}

layer_rows <- vector("list", length(count_layers))
for (index in seq_along(count_layers)) {
  layer <- count_layers[[index]]
  layer_cells <- Cells(object[["RNA"]], layer = layer)
  cell_index <- match(layer_cells, assignments$Cells)
  dataset <- unique(as.character(object$Dataset[cell_index]))
  if (anyNA(cell_index) || anyDuplicated(layer_cells) ||
    length(dataset) != 1L || is.na(dataset)) {
    stop("Lineage count layer does not map to one source dataset: ", layer)
  }
  layer_rows[[index]] <- data.table(
    Dataset = dataset, Layer = layer, Full_Cells = length(layer_cells),
    Eligible_Cells = sum(assignments$Subtype_Inference_Eligible[cell_index])
  )
}
layer_map <- rbindlist(layer_rows)
if (anyDuplicated(layer_map$Dataset) ||
  sum(layer_map$Full_Cells) != nrow(assignments) ||
  sum(layer_map$Eligible_Cells) != nrow(marker_assignments) ||
  any(layer_map$Eligible_Cells <= 0L)) {
  stop("Lineage dataset-layer coverage is incomplete")
}
setorder(layer_map, Dataset)

safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
checkpoint_files <- stats::setNames(
  file.path(
    checkpoint_dir,
    paste0(safe_name(layer_map$Dataset), "_aggregate.rds")
  ),
  layer_map$Dataset
)

for (layer_index in seq_len(nrow(layer_map))) {
  dataset <- layer_map$Dataset[[layer_index]]
  layer <- layer_map$Layer[[layer_index]]
  checkpoint_file <- checkpoint_files[[dataset]]
  if (file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    if (!identical(checkpoint$Dataset, dataset) ||
      !identical(checkpoint$Layer, layer) ||
      !identical(checkpoint$Cluster_Levels, cluster_levels)) {
      stop("Lineage marker checkpoint contract changed: ", dataset)
    }
    next
  }
  message("[lineage-markers] aggregating ", lineage_id, " dataset ", dataset)
  layer_cells <- Cells(object[["RNA"]], layer = layer)
  cell_index <- match(layer_cells, assignments$Cells)
  eligible_index <- assignments$Subtype_Inference_Eligible[cell_index]
  layer_cells <- layer_cells[eligible_index]
  cell_index <- cell_index[eligible_index]
  layer_clusters <- factor(
    assignments$Lineage_Cluster[cell_index],
    levels = cluster_levels
  )
  cluster_sizes <- tabulate(
    as.integer(layer_clusters),
    nbins = length(cluster_levels)
  )
  membership <- sparseMatrix(
    i = seq_along(layer_cells), j = as.integer(layer_clusters),
    x = 1, dims = c(length(layer_cells), length(cluster_levels)),
    dimnames = list(layer_cells, cluster_levels)
  )
  layer_features <- intersect(
    features, Features(object[["RNA"]], layer = layer)
  )
  counts <- LayerData(
    object,
    assay = "RNA", layer = layer,
    cells = layer_cells, features = layer_features, fast = FALSE
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  if (!identical(colnames(counts), layer_cells) ||
    any(!is.finite(counts@x)) || any(counts@x < 0)) {
    stop("Lineage count layer failed value or order checks: ", dataset)
  }
  library_sizes <- Matrix::colSums(counts)
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop("Lineage count layer contains a non-positive library: ", dataset)
  }
  normalized <- counts
  scale <- normalization_scale / library_sizes
  normalized@x <- log1p(
    normalized@x * rep.int(scale, diff(normalized@p))
  )
  sum_log_normalized <- as.matrix(normalized %*% membership)
  detected <- counts
  detected@x[] <- 1
  detected_cells <- as.matrix(detected %*% membership)
  checkpoint <- list(
    Dataset = dataset,
    Layer = layer,
    Cells = length(layer_cells),
    Features = rownames(sum_log_normalized),
    Cluster_Levels = cluster_levels,
    Cluster_Sizes = cluster_sizes,
    Sum_LogNormalized = sum_log_normalized,
    Detected_Cells = detected_cells
  )
  temporary <- paste0(checkpoint_file, ".tmp.", Sys.getpid())
  saveRDS(checkpoint, temporary, compress = TRUE)
  if (!file.rename(temporary, checkpoint_file)) {
    unlink(temporary)
    stop("Could not publish lineage marker checkpoint: ", dataset)
  }
  rm(
    counts, normalized, detected, membership, checkpoint,
    sum_log_normalized, detected_cells
  )
  gc()
}

n_features <- length(features)
n_clusters <- length(cluster_levels)
n_datasets <- nrow(layer_map)
dataset_effect <- array(
  NA_real_,
  dim = c(n_features, n_clusters, n_datasets),
  dimnames = list(features, cluster_levels, layer_map$Dataset)
)
dataset_pct_difference <- dataset_effect
global_sum <- matrix(
  0, n_features, n_clusters,
  dimnames = list(features, cluster_levels)
)
global_detected <- global_sum
global_cluster_cells <- global_sum
global_available_cells <- stats::setNames(numeric(n_features), features)
dataset_audit <- vector("list", n_datasets)

for (dataset_index in seq_len(n_datasets)) {
  dataset <- layer_map$Dataset[[dataset_index]]
  checkpoint <- readRDS(checkpoint_files[[dataset]])
  feature_index <- match(checkpoint$Features, features)
  cluster_sizes <- as.numeric(checkpoint$Cluster_Sizes)
  total_cells <- sum(cluster_sizes)
  valid_clusters <- cluster_sizes >= minimum_cluster_cells_per_dataset &
    (total_cells - cluster_sizes) >= minimum_cluster_cells_per_dataset
  sum_log <- checkpoint$Sum_LogNormalized
  detected <- checkpoint$Detected_Cells
  mean_in <- sweep(sum_log, 2L, pmax(cluster_sizes, 1), "/")
  mean_out <- sweep(
    rowSums(sum_log) - sum_log, 2L, total_cells - cluster_sizes, "/"
  )
  pct_in <- sweep(detected, 2L, pmax(cluster_sizes, 1), "/")
  pct_out <- sweep(
    rowSums(detected) - detected, 2L, total_cells - cluster_sizes, "/"
  )
  dataset_effect[feature_index, valid_clusters, dataset_index] <-
    (mean_in - mean_out)[, valid_clusters, drop = FALSE]
  dataset_pct_difference[feature_index, valid_clusters, dataset_index] <-
    (pct_in - pct_out)[, valid_clusters, drop = FALSE]
  global_sum[feature_index, ] <- global_sum[feature_index, ] + sum_log
  global_detected[feature_index, ] <- global_detected[feature_index, ] + detected
  global_cluster_cells[feature_index, ] <-
    global_cluster_cells[feature_index, ] + matrix(
      rep(cluster_sizes, each = length(feature_index)),
      nrow = length(feature_index), ncol = n_clusters
    )
  global_available_cells[feature_index] <-
    global_available_cells[feature_index] + total_cells
  dataset_audit[[dataset_index]] <- data.table(
    Dataset = dataset, Cells = total_cells,
    Features = length(feature_index),
    Clusters_With_At_Least_20_Cells = sum(valid_clusters)
  )
  rm(checkpoint, sum_log, detected, mean_in, mean_out, pct_in, pct_out)
  gc()
}

global_mean_in <- global_sum / global_cluster_cells
global_pct_in <- global_detected / global_cluster_cells
global_mean_out <- (rowSums(global_sum) - global_sum) /
  (global_available_cells - global_cluster_cells)
global_pct_out <- (rowSums(global_detected) - global_detected) /
  (global_available_cells - global_cluster_cells)
marker_rows <- vector("list", n_clusters)
top_rows <- vector("list", n_clusters * 2L)
top_index <- 0L
for (cluster_index in seq_len(n_clusters)) {
  effect <- matrix(
    dataset_effect[, cluster_index, ],
    nrow = n_features,
    ncol = n_datasets
  )
  pct_difference <- matrix(
    dataset_pct_difference[, cluster_index, ],
    nrow = n_features,
    ncol = n_datasets
  )
  evidence <- rowSums(is.finite(effect))
  mean_effect <- rowMeans2(effect, na.rm = TRUE)
  median_effect <- rowMedians(effect, na.rm = TRUE)
  mean_pct <- rowMeans2(pct_difference, na.rm = TRUE)
  positive_fraction <- rowSums(effect > 0, na.rm = TRUE) / evidence
  negative_fraction <- rowSums(effect < 0, na.rm = TRUE) / evidence
  zero <- evidence == 0L
  mean_effect[zero] <- median_effect[zero] <- mean_pct[zero] <- NA_real_
  positive_fraction[zero] <- negative_fraction[zero] <- NA_real_
  current <- data.table(
    Lineage_Cluster = cluster_levels[[cluster_index]],
    Gene = features,
    Global_Mean_LogNormalized_In = global_mean_in[, cluster_index],
    Global_Mean_LogNormalized_Out = global_mean_out[, cluster_index],
    Global_Pct_Expressed_In = global_pct_in[, cluster_index],
    Global_Pct_Expressed_Out = global_pct_out[, cluster_index],
    Evidence_Datasets = evidence,
    Dataset_Equal_Mean_LogNormalized_Difference = mean_effect,
    Dataset_Median_LogNormalized_Difference = median_effect,
    Dataset_Equal_Mean_Pct_Expressed_Difference = mean_pct,
    Positive_Dataset_Fraction = positive_fraction,
    Negative_Dataset_Fraction = negative_fraction
  )
  marker_rows[[cluster_index]] <- current
  positive <- head(current[
    Evidence_Datasets > 0 & Dataset_Equal_Mean_LogNormalized_Difference > 0
  ][order(
    -Dataset_Equal_Mean_LogNormalized_Difference,
    -Positive_Dataset_Fraction,
    -Dataset_Equal_Mean_Pct_Expressed_Difference,
    Gene
  )], top_markers_per_direction)
  negative <- head(current[
    Evidence_Datasets > 0 & Dataset_Equal_Mean_LogNormalized_Difference < 0
  ][order(
    Dataset_Equal_Mean_LogNormalized_Difference,
    -Negative_Dataset_Fraction,
    Dataset_Equal_Mean_Pct_Expressed_Difference,
    Gene
  )], top_markers_per_direction)
  positive[, `:=`(Direction = "positive", Rank = seq_len(.N))]
  negative[, `:=`(Direction = "negative", Rank = seq_len(.N))]
  top_index <- top_index + 1L
  top_rows[[top_index]] <- positive
  top_index <- top_index + 1L
  top_rows[[top_index]] <- negative
}
marker_statistics <- rbindlist(marker_rows)
top_markers <- rbindlist(top_rows, use.names = TRUE, fill = TRUE)
setcolorder(
  top_markers,
  c("Lineage_Cluster", "Direction", "Rank", "Gene", setdiff(
    names(top_markers), c("Lineage_Cluster", "Direction", "Rank", "Gene")
  ))
)

metadata_fields <- intersect(
  c(
    "Dataset", "Global_Donor_ID", "AgeIntervalID", "BrainRegion", "Sex",
    "Donor_ID", "Source_CellType", "Source_CellType_Original_Label",
    "source_cell_type_original_label", "source_cell_type_label",
    "source_cell_type_level_1", "source_cell_type_level_2",
    "source_cell_type_level_3", "source_cluster_id"
  ),
  names(object[[]])
)
metadata <- as.data.table(object[[]][, metadata_fields, drop = FALSE])
metadata[, Lineage_Cluster := assignments$Lineage_Cluster]
metadata[, Global_Cluster := assignments$Global_Cluster]
metadata[, Subtype_Inference_Eligible :=
  assignments$Subtype_Inference_Eligible]
composition <- rbindlist(lapply(metadata_fields, function(field) {
  values <- normalize_missing_metadata(metadata[[field]])
  values[is.na(values)] <- "<missing>"
  data.table(
    Lineage_Cluster = metadata$Lineage_Cluster,
    Field = field,
    Value = values
  )[, .(Cells = .N), by = .(Lineage_Cluster, Field, Value)][
    , Cluster_Fraction := Cells / sum(Cells),
    by = Lineage_Cluster
  ]
}))

eligibility_by_cluster <- assignments[, .(
  Full_Cells = .N,
  Subtype_Inference_Eligible_Cells = sum(Subtype_Inference_Eligible),
  Subtype_Inference_Excluded_Cells = sum(!Subtype_Inference_Eligible)
), by = .(Lineage_Cluster)][
  , Eligible_Fraction := Subtype_Inference_Eligible_Cells / Full_Cells
]
setorder(eligibility_by_cluster, Lineage_Cluster)
subtype_exclusions <- main_celltype_audit[
  Cluster %chin% excluded_subtype_clusters,
  .(Cluster, Main_CellType, Exclude_From_Subtype_Inference)
]
if (nrow(subtype_exclusions) > 0L) {
  excluded_counts <- assignments[
    Subtype_Inference_Eligible == FALSE,
    .(Excluded_Cells = .N),
    by = .(Cluster = Global_Cluster)
  ]
  subtype_exclusions <- merge(
    subtype_exclusions, excluded_counts,
    by = "Cluster", all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(subtype_exclusions$Excluded_Cells)) {
    stop("Subtype-inference exclusion counts are incomplete")
  }
} else {
  subtype_exclusions[, Excluded_Cells := integer()]
}

statistics_file <- file.path(output_dir, "lineage_marker_statistics.tsv.gz")
top_file <- file.path(output_dir, "lineage_top_positive_negative_markers.tsv")
composition_file <- file.path(output_dir, "lineage_cluster_composition.tsv.gz")
dataset_audit_file <- file.path(output_dir, "lineage_dataset_marker_audit.tsv")
eligibility_file <- file.path(
  output_dir, "lineage_cluster_subtype_eligibility.tsv"
)
exclusion_file <- file.path(
  output_dir, "lineage_subtype_inference_exclusions.tsv"
)
write_tsv(as.data.frame(marker_statistics), statistics_file)
write_tsv(as.data.frame(top_markers), top_file)
write_tsv(as.data.frame(composition), composition_file)
write_tsv(as.data.frame(rbindlist(dataset_audit)), dataset_audit_file)
write_tsv(as.data.frame(eligibility_by_cluster), eligibility_file)
write_tsv(as.data.frame(subtype_exclusions), exclusion_file)
output_files <- c(
  statistics_file, top_file, composition_file, dataset_audit_file,
  eligibility_file, exclusion_file
)
write_tsv(
  data.frame(
    File = output_files,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "lineage_marker_output_manifest.tsv")
)
write_tsv(
  data.frame(
    Field = c(
      "Lineage_ID", "Full_Cells", "Subtype_Inference_Eligible_Cells",
      "Subtype_Inference_Excluded_Cells", "Subtype_Eligible_Clusters",
      "Genes", "Datasets", "Normalization",
      "Minimum_Cluster_Cells_Per_Dataset", "Dataset_Weighting",
      "Cell_Level_Hypothesis_Test", "Positive_And_Negative_Markers"
    ),
    Value = c(
      lineage_id, nrow(assignments), nrow(marker_assignments),
      nrow(assignments) - nrow(marker_assignments), n_clusters,
      n_features, n_datasets,
      "log1p(10000 * raw count / cell library size)",
      minimum_cluster_cells_per_dataset,
      "source datasets contribute equally to marker-effect summaries",
      "FALSE; cell-level expression fractions are descriptive",
      "full statistics plus top 50 per direction per lineage cluster"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "lineage_marker_contract.tsv")
)
message(
  "[lineage-markers] completed ", lineage_id, " markers for ",
  n_clusters, " subtype-eligible clusters and ",
  nrow(marker_assignments), " eligible cells (",
  nrow(assignments), " full lineage cells retained in clustering)"
)
