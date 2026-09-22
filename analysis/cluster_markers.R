suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(matrixStats)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")


integration_dir <- brainomics_data_path("integration_25")
annotation_dir <- file.path(integration_dir, "annotation")
source("functions/utils.R")
assignment_column <- "Cluster"
output_dir <- file.path(annotation_dir, "markers")
checkpoint_dir <- file.path(output_dir, "dataset_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
object_file <- file.path(integration_dir, "objects_merged.rds")
feature_file <- file.path(annotation_dir, "rna_features.tsv.gz")
normalization_scale <- 10000
minimum_cluster_cells_per_dataset <- 20L
top_markers_per_direction <- 25L
assignments <- as.data.table(read_celltype_assignments())
assignments[, Cluster := as.character(as.integer(sub("^C", "", Cluster)))]
cluster_levels <- as.character(0:74)

feature_table <- read.delim(
  gzfile(feature_file),
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!identical(names(feature_table), c("Feature_Order", "Feature")) ||
  !identical(feature_table$Feature_Order, seq_len(nrow(feature_table))) ||
  anyNA(feature_table$Feature) ||
  anyDuplicated(feature_table$Feature)) {
  stop("Final RNA feature table is invalid")
}
features <- as.character(feature_table$Feature)

thisutils::log_message("[cluster-markers] ", "Loading the complete merged count-layer object")
objects <- readRDS(object_file)
if (!inherits(objects, "Seurat") ||
  !identical(colnames(objects), assignments$Cells)) {
  stop("Merged count-layer object differs from the clustered cell cohort")
}
if (!"Dataset" %in% names(objects[[]])) {
  stop("Merged count-layer object lacks Dataset metadata")
}

count_layers <- Layers(objects, assay = "RNA", search = "^counts")
if (length(count_layers) != length(reference_datasets())) {
  stop("Merged count-layer object does not contain one count layer per dataset")
}

layer_rows <- vector("list", length(count_layers))
for (layer_index in seq_along(count_layers)) {
  layer <- count_layers[[layer_index]]
  layer_cells <- Cells(objects[["RNA"]], layer = layer)
  cell_index <- match(layer_cells, assignments$Cells)
  if (anyNA(cell_index) || anyDuplicated(layer_cells)) {
    stop(layer, " contains unmatched or duplicated cells")
  }
  datasets <- unique(as.character(objects$Dataset[cell_index]))
  if (length(datasets) != 1L || is.na(datasets)) {
    stop(layer, " does not map to exactly one source dataset")
  }
  layer_rows[[layer_index]] <- data.frame(
    Dataset = datasets,
    Layer = layer,
    Cells = length(layer_cells),
    stringsAsFactors = FALSE
  )
}
layer_map <- rbindlist(layer_rows)
if (anyDuplicated(layer_map$Dataset) ||
  !setequal(layer_map$Dataset, reference_datasets()) ||
  sum(layer_map$Cells) != nrow(assignments)) {
  stop("Count-layer to dataset mapping differs from the reference contract")
}
layer_map[, Dataset_Order := match(Dataset, reference_datasets())]
setorder(layer_map, Dataset_Order)
layer_map[, Dataset_Order := NULL]

safe_name <- function(value) {
  gsub("[^A-Za-z0-9_.-]", "_", value)
}

checkpoint_is_valid <- function(checkpoint, dataset, layer, cells) {
  is.list(checkpoint) &&
    identical(checkpoint$Dataset, dataset) &&
    identical(checkpoint$Layer, layer) &&
    identical(as.numeric(checkpoint$Cells), as.numeric(cells)) &&
    identical(checkpoint$Cluster_Levels, cluster_levels) &&
    identical(
      as.numeric(checkpoint$Normalization_Scale),
      normalization_scale
    ) &&
    is.matrix(checkpoint$Sum_LogNormalized) &&
    is.matrix(checkpoint$Detected_Cells) &&
    ncol(checkpoint$Sum_LogNormalized) == length(cluster_levels) &&
    identical(
      dim(checkpoint$Sum_LogNormalized),
      dim(checkpoint$Detected_Cells)
    ) &&
    length(checkpoint$Features) == nrow(checkpoint$Sum_LogNormalized) &&
    length(checkpoint$Cluster_Sizes) == length(cluster_levels)
}

checkpoint_files <- stats::setNames(
  file.path(
    checkpoint_dir,
    paste0(safe_name(layer_map$Dataset), "_marker_aggregate.rds")
  ),
  layer_map$Dataset
)

for (layer_index in seq_len(nrow(layer_map))) {
  dataset <- layer_map$Dataset[[layer_index]]
  layer <- layer_map$Layer[[layer_index]]
  checkpoint_file <- checkpoint_files[[dataset]]
  if (file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    if (checkpoint_is_valid(
      checkpoint,
      dataset,
      layer,
      layer_map$Cells[[layer_index]]
    )) {
      thisutils::log_message("[cluster-markers] ", paste("Using completed marker aggregate for", dataset))
      next
    }
    stop("Existing marker aggregate differs from its contract: ", dataset)
  }

  thisutils::log_message("[cluster-markers] ", paste("Aggregating complete count layer for", dataset))
  layer_cells <- Cells(objects[["RNA"]], layer = layer)
  cell_index <- match(layer_cells, assignments$Cells)
  layer_clusters <- factor(
    as.character(assignments$Cluster[cell_index]),
    levels = cluster_levels
  )
  if (anyNA(layer_clusters)) {
    stop(dataset, " has cells without a primary cluster assignment")
  }
  cluster_sizes <- tabulate(
    as.integer(layer_clusters),
    nbins = length(cluster_levels)
  )
  membership <- sparseMatrix(
    i = seq_along(layer_cells),
    j = as.integer(layer_clusters),
    x = rep.int(1, length(layer_cells)),
    dims = c(length(layer_cells), length(cluster_levels)),
    dimnames = list(layer_cells, cluster_levels)
  )

  layer_features <- intersect(
    features,
    Features(objects[["RNA"]], layer = layer)
  )
  counts <- LayerData(
    objects,
    assay = "RNA",
    layer = layer,
    cells = layer_cells,
    features = layer_features,
    fast = FALSE
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  if (!identical(colnames(counts), layer_cells) ||
    anyDuplicated(rownames(counts)) ||
    any(!is.finite(counts@x)) ||
    any(counts@x < 0)) {
    stop(dataset, " count layer failed matrix or order checks")
  }
  library_sizes <- Matrix::colSums(counts)
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop(dataset, " contains a cell with a non-positive RNA library size")
  }

  normalized <- counts
  column_scale <- normalization_scale / library_sizes
  normalized@x <- log1p(
    normalized@x * rep.int(column_scale, diff(normalized@p))
  )
  sum_log_normalized <- as.matrix(normalized %*% membership)
  rm(normalized, column_scale)
  gc()

  detected <- counts
  detected@x[] <- 1
  detected_cells <- as.matrix(detected %*% membership)
  rm(detected, counts, membership)
  gc()

  checkpoint <- list(
    Dataset = dataset,
    Layer = layer,
    Cells = length(layer_cells),
    Features = rownames(sum_log_normalized),
    Cluster_Levels = cluster_levels,
    Cluster_Sizes = cluster_sizes,
    Normalization = "log1p(10000 * raw count / cell library size)",
    Normalization_Scale = normalization_scale,
    Sum_LogNormalized = sum_log_normalized,
    Detected_Cells = detected_cells
  )
  temporary_checkpoint <- paste0(checkpoint_file, ".tmp.", Sys.getpid())
  saveRDS(checkpoint, temporary_checkpoint, compress = TRUE)
  if (!file.rename(temporary_checkpoint, checkpoint_file)) {
    unlink(temporary_checkpoint)
    stop("Could not publish marker aggregate for ", dataset)
  }
  rm(checkpoint, sum_log_normalized, detected_cells)
  gc()
}

thisutils::log_message("[cluster-markers] ", "Combining dataset-level marker aggregates")
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
  0,
  nrow = n_features,
  ncol = n_clusters,
  dimnames = list(features, cluster_levels)
)
global_detected <- global_sum
global_cluster_cells <- global_sum
global_available_cells <- stats::setNames(numeric(n_features), features)
dataset_audit_rows <- vector("list", n_datasets)

for (dataset_index in seq_len(n_datasets)) {
  dataset <- layer_map$Dataset[[dataset_index]]
  checkpoint <- readRDS(checkpoint_files[[dataset]])
  if (!checkpoint_is_valid(
    checkpoint,
    dataset,
    layer_map$Layer[[dataset_index]],
    layer_map$Cells[[dataset_index]]
  )) {
    stop("Marker aggregate failed reload validation: ", dataset)
  }
  feature_index <- match(checkpoint$Features, features)
  if (anyNA(feature_index) || anyDuplicated(feature_index)) {
    stop(dataset, " marker aggregate contains invalid features")
  }
  cluster_sizes <- as.numeric(checkpoint$Cluster_Sizes)
  total_cells <- sum(cluster_sizes)
  valid_clusters <- cluster_sizes >= minimum_cluster_cells_per_dataset &
    (total_cells - cluster_sizes) >= minimum_cluster_cells_per_dataset

  sum_log <- checkpoint$Sum_LogNormalized
  detected <- checkpoint$Detected_Cells
  total_sum_log <- rowSums(sum_log)
  total_detected <- rowSums(detected)
  mean_in <- sweep(sum_log, 2L, pmax(cluster_sizes, 1), "/")
  mean_out <- sweep(
    total_sum_log - sum_log,
    2L,
    total_cells - cluster_sizes,
    "/"
  )
  pct_in <- sweep(detected, 2L, pmax(cluster_sizes, 1), "/")
  pct_out <- sweep(
    total_detected - detected,
    2L,
    total_cells - cluster_sizes,
    "/"
  )
  dataset_effect[feature_index, valid_clusters, dataset_index] <-
    (mean_in - mean_out)[, valid_clusters, drop = FALSE]
  dataset_pct_difference[feature_index, valid_clusters, dataset_index] <-
    (pct_in - pct_out)[, valid_clusters, drop = FALSE]

  global_sum[feature_index, ] <- global_sum[feature_index, ] + sum_log
  global_detected[feature_index, ] <-
    global_detected[feature_index, ] + detected
  global_cluster_cells[feature_index, ] <-
    global_cluster_cells[feature_index, ] + matrix(
      rep(cluster_sizes, each = length(feature_index)),
      nrow = length(feature_index),
      ncol = n_clusters
    )
  global_available_cells[feature_index] <-
    global_available_cells[feature_index] + total_cells

  dataset_audit_rows[[dataset_index]] <- data.frame(
    Dataset = dataset,
    Layer = checkpoint$Layer,
    Cells = total_cells,
    Features = length(checkpoint$Features),
    Nonzero_Cluster_Feature_Values = sum(detected > 0),
    Clusters_With_At_Least_20_Cells = sum(valid_clusters),
    Checkpoint = checkpoint_files[[dataset]],
    stringsAsFactors = FALSE
  )
  rm(
    checkpoint,
    sum_log,
    detected,
    mean_in,
    mean_out,
    pct_in,
    pct_out,
    total_sum_log,
    total_detected
  )
  gc()
}

global_mean_in <- global_sum / global_cluster_cells
global_pct_in <- global_detected / global_cluster_cells
global_sum_all_clusters <- rowSums(global_sum)
global_detected_all_clusters <- rowSums(global_detected)
global_mean_out <- (global_sum_all_clusters - global_sum) /
  (global_available_cells - global_cluster_cells)
global_pct_out <- (global_detected_all_clusters - global_detected) /
  (global_available_cells - global_cluster_cells)

marker_rows <- vector("list", n_clusters)
top_rows <- vector("list", n_clusters * 2L)
top_index <- 0L
for (cluster_index in seq_len(n_clusters)) {
  cluster <- cluster_levels[[cluster_index]]
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
  finite_effect <- is.finite(effect)
  evidence_datasets <- rowSums(finite_effect)
  mean_effect <- rowMeans2(effect, na.rm = TRUE)
  median_effect <- rowMedians(effect, na.rm = TRUE)
  iqr_effect <- rowIQRs(effect, na.rm = TRUE)
  mean_pct_difference <- rowMeans2(pct_difference, na.rm = TRUE)
  positive_fraction <- rowSums(effect > 0, na.rm = TRUE) / evidence_datasets
  negative_fraction <- rowSums(effect < 0, na.rm = TRUE) / evidence_datasets
  zero_evidence <- evidence_datasets == 0L
  mean_effect[zero_evidence] <- NA_real_
  median_effect[zero_evidence] <- NA_real_
  iqr_effect[zero_evidence] <- NA_real_
  mean_pct_difference[zero_evidence] <- NA_real_
  positive_fraction[zero_evidence] <- NA_real_
  negative_fraction[zero_evidence] <- NA_real_

  marker_rows[[cluster_index]] <- data.table(
    Cluster = cluster,
    Gene = features,
    Global_Mean_LogNormalized_In = global_mean_in[, cluster_index],
    Global_Mean_LogNormalized_Out = global_mean_out[, cluster_index],
    Global_LogNormalized_Difference =
      global_mean_in[, cluster_index] - global_mean_out[, cluster_index],
    Global_Pct_Expressed_In = global_pct_in[, cluster_index],
    Global_Pct_Expressed_Out = global_pct_out[, cluster_index],
    Global_Pct_Expressed_Difference =
      global_pct_in[, cluster_index] - global_pct_out[, cluster_index],
    Evidence_Datasets = evidence_datasets,
    Dataset_Equal_Mean_LogNormalized_Difference = mean_effect,
    Dataset_Median_LogNormalized_Difference = median_effect,
    Dataset_IQR_LogNormalized_Difference = iqr_effect,
    Dataset_Equal_Mean_Pct_Expressed_Difference = mean_pct_difference,
    Positive_Dataset_Fraction = positive_fraction,
    Negative_Dataset_Fraction = negative_fraction
  )
  current <- marker_rows[[cluster_index]][Evidence_Datasets > 0L]
  positive <- head(current[
    Dataset_Equal_Mean_LogNormalized_Difference > 0
  ][order(
    -Dataset_Equal_Mean_LogNormalized_Difference,
    -Positive_Dataset_Fraction,
    -Global_Pct_Expressed_Difference,
    Gene
  )], top_markers_per_direction)
  negative <- head(current[
    Dataset_Equal_Mean_LogNormalized_Difference < 0
  ][order(
    Dataset_Equal_Mean_LogNormalized_Difference,
    -Negative_Dataset_Fraction,
    Global_Pct_Expressed_Difference,
    Gene
  )], top_markers_per_direction)
  positive[, `:=`(
    Direction = "positive",
    Rank = seq_len(.N)
  )]
  negative[, `:=`(
    Direction = "negative",
    Rank = seq_len(.N)
  )]
  top_index <- top_index + 1L
  top_rows[[top_index]] <- positive
  top_index <- top_index + 1L
  top_rows[[top_index]] <- negative
}
marker_statistics <- rbindlist(marker_rows, use.names = TRUE)
top_markers <- rbindlist(top_rows, use.names = TRUE, fill = TRUE)
setcolorder(
  top_markers,
  c("Cluster", "Direction", "Rank", "Gene", setdiff(
    names(top_markers),
    c("Cluster", "Direction", "Rank", "Gene")
  ))
)

cluster_annotation <- as.data.table(read_tsv("results/annotation/cluster_annotation.tsv"))
cluster_annotation[, Cluster := as.character(as.integer(sub("^C", "", Cluster)))]
positive_strings <- top_markers[
  Direction == "positive" & Rank <= 20L,
  .(Positive_Markers = paste(Gene, collapse = ";")),
  by = Cluster
]
negative_strings <- top_markers[
  Direction == "negative" & Rank <= 20L,
  .(Negative_Markers = paste(Gene, collapse = ";")),
  by = Cluster
]
cluster_summary <- merge(
  cluster_annotation,
  positive_strings,
  by = "Cluster",
  all.x = TRUE,
  sort = FALSE
)
cluster_summary <- merge(
  cluster_summary,
  negative_strings,
  by = "Cluster",
  all.x = TRUE,
  sort = FALSE
)
cluster_summary[, Cluster_Order := as.integer(Cluster)]
setorder(cluster_summary, Cluster_Order)
cluster_summary[, Cluster_Order := NULL]

metadata_fields <- intersect(
  c("Dataset", "AgeIntervalID", "BrainRegion", "Sex", "Technology", "Modality"),
  names(objects[[]])
)
metadata <- objects[[]][, metadata_fields, drop = FALSE]
if (!identical(rownames(metadata), assignments$Cells)) {
  stop("Metadata and cluster assignment orders differ")
}
composition_rows <- vector("list", length(metadata_fields))
for (field_index in seq_along(metadata_fields)) {
  field <- metadata_fields[[field_index]]
  value <- normalize_missing_metadata(metadata[[field]])
  value[is.na(value)] <- "<missing>"
  composition <- data.table(
    Cluster = assignments$Cluster,
    Field = field,
    Value = value
  )[, .(Cells = .N), by = .(Cluster, Field, Value)]
  composition[, Cluster_Fraction := Cells / sum(Cells), by = Cluster]
  composition_rows[[field_index]] <- composition
}
metadata_composition <- rbindlist(composition_rows, use.names = TRUE)

statistics_file <- file.path(output_dir, "cluster_marker_statistics.tsv.gz")
top_file <- file.path(output_dir, "cluster_top_positive_negative_markers.tsv")
summary_file <- file.path(output_dir, "cluster_marker_summary.tsv")
composition_file <- file.path(output_dir, "cluster_metadata_composition.tsv.gz")
dataset_audit_file <- file.path(output_dir, "dataset_marker_audit.tsv")
write_tsv(as.data.frame(marker_statistics), statistics_file)
write_tsv(as.data.frame(top_markers), top_file)
write_tsv(as.data.frame(cluster_summary), summary_file)
write_tsv(as.data.frame(metadata_composition), composition_file)
write_tsv(as.data.frame(rbindlist(dataset_audit_rows)), dataset_audit_file)

output_manifest <- data.frame(
  File = basename(c(
    statistics_file,
    top_file,
    summary_file,
    composition_file,
    dataset_audit_file
  )),
  stringsAsFactors = FALSE
)
write_tsv(output_manifest, file.path(output_dir, "marker_output_manifest.tsv"))
write_tsv(
  data.frame(
    Field = c(
      "Cells", "Clusters", "Genes", "Datasets", "Input_Object",
      "Input_Object_Size_Bytes",
      "Normalization", "Minimum_Cluster_Cells_Per_Dataset",
      "Dataset_Weighting", "Cell_Level_Hypothesis_Test",
      "Positive_And_Negative_Markers", "Assignment_Column",
      "R_Version", "Matrix_Version", "SeuratObject_Version",
      "matrixStats_Version"
    ),
    Value = c(
      nrow(assignments), n_clusters, n_features, n_datasets, basename(object_file),
      file.info(object_file)$size,
      "log1p(10000 * raw count / cell library size)",
      minimum_cluster_cells_per_dataset,
      "source datasets contribute equally to marker effect summaries",
      "FALSE; cell-level expression fractions are descriptive",
      "full statistics plus top 25 per direction per cluster",
      assignment_column, R.version.string,
      as.character(packageVersion("Matrix")),
      as.character(packageVersion("SeuratObject")),
      as.character(packageVersion("matrixStats"))
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "marker_contract.tsv")
)

thisutils::log_message("[cluster-markers] ", paste(
  "Completed positive and negative marker evidence for",
  n_clusters,
  "clusters across",
  n_datasets,
  "datasets and",
  format(n_features, big.mark = ","),
  "genes"
))
