suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/metadata_schema.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")


res_dir <- brainomics_data_path("integration_25")
evaluation_dir <- file.path(res_dir, "evaluation")
checkpoint_dir <- file.path(
  evaluation_dir,
  "lisi_checkpoints_double_precision"
)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

object_file <- file.path(res_dir, "objects_integrated.rds")
output_file <- file.path(res_dir, "lisi_results.rds")
overwrite <- tolower(Sys.getenv(
  "BRAINOMICS_EVALUATION_OVERWRITE",
  unset = "false"
)) %in% c("true", "t", "1")
requested_method <- trimws(Sys.getenv(
  "BRAINOMICS_EVALUATION_METHOD",
  unset = "all"
))
perplexity <- as.numeric(Sys.getenv(
  "BRAINOMICS_LISI_PERPLEXITY",
  unset = "30"
))
nn_eps <- as.numeric(Sys.getenv(
  "BRAINOMICS_LISI_NN_EPS",
  unset = "0.1"
))

if (!file.exists(object_file)) {
  stop("Final integrated reference is missing: ", object_file)
}
if (!requireNamespace("lisi", quietly = TRUE)) {
  stop("The pinned lisi package is required for integration evaluation")
}
lisi_contract <- read_lisi_numeric_contract()

reduction_contract_file <- file.path(
  res_dir,
  "integration_reduction_audit.tsv"
)
reduction_contract <- read_integration_reduction_cell_contract(
  reduction_contract_file
)
cell_hash <- reduction_contract$Cell_Order_SHA256

thisutils::log_message("[integration-05] ", "Loading the final four-method reference")
object <- load_processed_object(object_file)
meta <- object[[]]
brain_region_before <- normalize_missing_metadata(meta$BrainRegion)
meta$BrainRegion <- standardize_brain_region(brain_region_before)
brain_region_changed <- !is.na(brain_region_before) &
  !is.na(meta$BrainRegion) &
  brain_region_before != meta$BrainRegion
brain_region_case_audit <- if (any(brain_region_changed)) {
  aggregate(
    rep(1L, sum(brain_region_changed)),
    by = list(
      Dataset = as.character(meta$Dataset[brain_region_changed]),
      Source_Label = brain_region_before[brain_region_changed],
      Harmonized_Label = meta$BrainRegion[brain_region_changed]
    ),
    FUN = sum
  )
} else {
  data.frame(
    Dataset = character(), Source_Label = character(),
    Harmonized_Label = character(), x = integer(),
    stringsAsFactors = FALSE
  )
}
names(brain_region_case_audit)[names(brain_region_case_audit) == "x"] <-
  "Cells"
write_tsv(
  brain_region_case_audit,
  file.path(evaluation_dir, "brain_region_case_canonicalization_audit.tsv")
)
method_identity <- integration_method_identity()
latent_identity <- method_identity[
  method_identity$Space == "latent", ,
  drop = FALSE
]
latent_identity <- latent_identity[match(
  method_levels,
  latent_identity$Display_Name
), , drop = FALSE]
valid_methods <- as.character(latent_identity$Display_Name)
if (!requested_method %in% c("all", "aggregate", valid_methods)) {
  stop(
    "Unknown BRAINOMICS_EVALUATION_METHOD: ",
    requested_method,
    "; expected all, aggregate, or one of ",
    paste(valid_methods, collapse = ", ")
  )
}
methods_to_compute <- if (requested_method %in% valid_methods) {
  requested_method
} else if (identical(requested_method, "aggregate")) {
  character()
} else {
  valid_methods
}

required_metadata <- c(
  "Dataset", "Global_Donor_ID", "AgeIntervalID", "BrainRegion",
  "Assay_Type", "Sequencing_Platform", "Library_Chemistry", "Sex"
)
missing_metadata <- setdiff(required_metadata, names(meta))
if (length(missing_metadata) > 0L) {
  stop(
    "Integration evaluation metadata is missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}
if (!identical(rownames(meta), colnames(object))) {
  stop("Final object metadata differs from matrix cell order")
}
if (nrow(meta) != reduction_contract$Cells) {
  stop("Final object cell count differs from integration reduction audit")
}
if (anyNA(meta$Dataset) || any(meta$Dataset == "")) {
  stop("Dataset labels must be complete for iLISI")
}

score_list <- vector("list", length(valid_methods))
names(score_list) <- valid_methods
for (display_name in methods_to_compute) {
  index <- match(display_name, latent_identity$Display_Name)
  reduction <- latent_identity$Reduction[[index]]
  checkpoint <- file.path(
    checkpoint_dir,
    paste0(tolower(display_name), "_dataset_lisi.rds")
  )
  if (file.exists(checkpoint) && !overwrite) {
    thisutils::log_message("[integration-05] ", paste("Resuming", display_name, "iLISI checkpoint"))
    value <- readRDS(checkpoint)
  } else {
    thisutils::log_message("[integration-05] ", paste("Computing full-cell dataset iLISI for", display_name))
    embedding <- Embeddings(object, reduction = reduction)
    if (!identical(rownames(embedding), rownames(meta)) ||
      any(!is.finite(embedding))) {
      stop(display_name, " latent embedding failed order or finite checks")
    }
    value <- lisi::compute_lisi(
      X = embedding,
      meta_data = meta[, "Dataset", drop = FALSE],
      label_colnames = "Dataset",
      perplexity = perplexity,
      nn_eps = nn_eps
    )$Dataset
    names(value) <- rownames(meta)
    value <- canonicalize_lisi_lower_bound(value)
    names(value) <- rownames(meta)
    temporary_checkpoint <- paste0(
      checkpoint,
      ".tmp.",
      Sys.getpid()
    )
    saveRDS(value, temporary_checkpoint, compress = TRUE)
    if (!file.rename(temporary_checkpoint, checkpoint)) {
      unlink(temporary_checkpoint)
      stop("Could not atomically publish ", display_name, " iLISI checkpoint")
    }
    rm(embedding)
    gc()
  }
  if (length(value) != nrow(meta) ||
    !identical(names(value), rownames(meta)) ||
    any(!is.finite(value)) ||
    any(value < 1)) {
    stop(display_name, " iLISI checkpoint failed validation")
  }
  score_list[[display_name]] <- value
  rm(value)
  gc()
}

if (requested_method %in% valid_methods) {
  thisutils::log_message("[integration-05] ", paste(
    "Completed full-cell latent-space iLISI checkpoint for",
    requested_method
  ))
  quit(save = "no", status = 0L)
}

for (display_name in valid_methods) {
  checkpoint <- file.path(
    checkpoint_dir,
    paste0(tolower(display_name), "_dataset_lisi.rds")
  )
  if (!file.exists(checkpoint)) {
    stop("Required iLISI checkpoint is missing: ", checkpoint)
  }
  value <- readRDS(checkpoint)
  value <- canonicalize_lisi_lower_bound(value)
  names(value) <- rownames(meta)
  if (length(value) != nrow(meta) ||
    !identical(names(value), rownames(meta)) ||
    any(!is.finite(value)) ||
    any(value < 1)) {
    stop(display_name, " iLISI checkpoint failed aggregation validation")
  }
  score_list[[display_name]] <- value
}

numeric_audit <- do.call(rbind, lapply(valid_methods, function(display_name) {
  value <- score_list[[display_name]]
  data.frame(
    Method = display_name,
    Numeric_Precision = "double",
    Pre_Adjustment_Minimum = attr(
      value,
      "pre_adjustment_minimum",
      exact = TRUE
    ),
    Theoretical_Lower_Bound_Adjustments = attr(
      value,
      "theoretical_lower_bound_adjustments",
      exact = TRUE
    ),
    Lower_Bound_Tolerance = attr(
      value,
      "lower_bound_tolerance",
      exact = TRUE
    ),
    stringsAsFactors = FALSE
  )
}))
write_tsv(numeric_audit, file.path(
  evaluation_dir,
  "lisi_numeric_audit.tsv"
))

lisi_results <- as.data.frame(score_list, check.names = FALSE)
rownames(lisi_results) <- rownames(meta)
attr(lisi_results, "method_identity") <- latent_identity
attr(lisi_results, "label_column") <- "Dataset"
attr(lisi_results, "space") <- "latent"
attr(lisi_results, "perplexity") <- perplexity
attr(lisi_results, "nn_eps") <- nn_eps
attr(lisi_results, "cell_order_sha256") <- cell_hash
attr(lisi_results, "numeric_precision") <- "double"
attr(lisi_results, "lisi_source_commit") <-
  unname(lisi_contract[["source_commit"]])
attr(lisi_results, "lisi_numeric_patch_sha256") <-
  unname(lisi_contract[["numeric_patch_sha256"]])
saveRDS(lisi_results, output_file, compress = TRUE)

group_means <- function(scores, group, group_data) {
  group <- normalize_missing_metadata(group)
  keep <- !is.na(group)
  if (!any(keep)) {
    return(NULL)
  }
  group <- factor(group[keep], levels = unique(group[keep]))
  sums <- rowsum(
    as.matrix(scores[keep, , drop = FALSE]),
    group = group,
    reorder = FALSE
  )
  cells <- as.numeric(rowsum(
    matrix(1, nrow = sum(keep), ncol = 1L),
    group = group,
    reorder = FALSE
  )[, 1L])
  means <- sweep(sums, 1L, cells, "/")
  first_index <- match(rownames(means), as.character(group))
  result <- cbind(
    group_data[which(keep)[first_index], , drop = FALSE],
    Cells = cells,
    as.data.frame(means, check.names = FALSE)
  )
  rownames(result) <- NULL
  result
}

known_donor <- normalize_missing_metadata(meta$Global_Donor_ID)
donor_group <- ifelse(
  is.na(known_donor),
  NA_character_,
  paste(meta$Dataset, known_donor, sep = "\r")
)
donor_scores <- group_means(
  lisi_results,
  donor_group,
  data.frame(
    Dataset = as.character(meta$Dataset),
    Global_Donor_ID = known_donor,
    stringsAsFactors = FALSE
  )
)
if (is.null(donor_scores)) {
  stop("No known donors are available for donor-level iLISI summaries")
}

dataset_scores <- group_means(
  donor_scores[, latent_identity$Display_Name, drop = FALSE],
  donor_scores$Dataset,
  data.frame(
    Dataset = donor_scores$Dataset,
    stringsAsFactors = FALSE
  )
)
names(dataset_scores)[names(dataset_scores) == "Cells"] <- "Known_Donors"

effect_rows <- lapply(latent_identity$Display_Name, function(method) {
  difference <- dataset_scores[[method]] - dataset_scores$Raw
  n <- length(difference)
  estimate <- mean(difference)
  standard_error <- if (n > 1L) stats::sd(difference) / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(0.975, df = n - 1L) else NA_real_
  data.frame(
    Method = method,
    Reference_Method = "Raw",
    Independent_Unit = "source dataset",
    Datasets = n,
    Mean_Paired_Difference = estimate,
    CI95_Lower = estimate - critical * standard_error,
    CI95_Upper = estimate + critical * standard_error,
    Median_Paired_Difference = stats::median(difference),
    IQR_Paired_Difference = stats::IQR(difference),
    stringsAsFactors = FALSE
  )
})
effects <- do.call(rbind, effect_rows)

stratified_rows <- list()
stratification_variables <- c(
  "AgeIntervalID", "BrainRegion", "Assay_Type",
  "Sequencing_Platform", "Library_Chemistry", "Sex"
)
for (variable in stratification_variables) {
  value <- normalize_missing_metadata(meta[[variable]])
  donor_stratum <- ifelse(
    is.na(known_donor) | is.na(value),
    NA_character_,
    paste(meta$Dataset, known_donor, value, sep = "\r")
  )
  donor_stratum_scores <- group_means(
    lisi_results,
    donor_stratum,
    data.frame(
      Dataset = as.character(meta$Dataset),
      Global_Donor_ID = known_donor,
      Stratum = value,
      stringsAsFactors = FALSE
    )
  )
  if (is.null(donor_stratum_scores)) {
    next
  }
  for (method in latent_identity$Display_Name) {
    split_values <- split(
      seq_len(nrow(donor_stratum_scores)),
      donor_stratum_scores$Stratum
    )
    stratified_rows[[paste(variable, method, sep = "\r")]] <- do.call(
      rbind,
      lapply(split_values, function(rows) {
        data.frame(
          Variable = variable,
          Stratum = donor_stratum_scores$Stratum[rows[[1L]]],
          Method = method,
          Cells = sum(donor_stratum_scores$Cells[rows]),
          Known_Donors = length(rows),
          Datasets = length(unique(donor_stratum_scores$Dataset[rows])),
          Donor_Equal_Mean_LISI = mean(
            donor_stratum_scores[[method]][rows]
          ),
          Donor_Equal_Median_LISI = stats::median(
            donor_stratum_scores[[method]][rows]
          ),
          stringsAsFactors = FALSE
        )
      })
    )
  }
}
stratified_scores <- do.call(rbind, stratified_rows)
rownames(stratified_scores) <- NULL

umap_identity <- method_identity[
  method_identity$Space == "UMAP", ,
  drop = FALSE
]
umap_identity <- umap_identity[match(
  latent_identity$Display_Name,
  umap_identity$Display_Name
), , drop = FALSE]
umap_plot_data <- data.frame(
  Cell = rownames(meta),
  Dataset = factor(meta$Dataset, levels = reference_datasets()),
  stringsAsFactors = FALSE
)
for (index in seq_len(nrow(umap_identity))) {
  method <- umap_identity$Display_Name[[index]]
  embedding <- Embeddings(
    object,
    reduction = umap_identity$Reduction[[index]]
  )
  if (!identical(rownames(embedding), rownames(meta)) ||
    ncol(embedding) != 2L || any(!is.finite(embedding))) {
    stop(method, " UMAP failed plotting-sidecar validation")
  }
  umap_plot_data[[paste0(method, "_1")]] <- embedding[, 1L]
  umap_plot_data[[paste0(method, "_2")]] <- embedding[, 2L]
}

write_tsv(donor_scores, file.path(
  evaluation_dir,
  "dataset_lisi_by_donor.tsv.gz"
))
write_tsv(dataset_scores, file.path(
  evaluation_dir,
  "dataset_lisi_dataset_equal.tsv"
))
write_tsv(effects, file.path(
  evaluation_dir,
  "dataset_lisi_effects.tsv"
))
write_tsv(stratified_scores, file.path(
  evaluation_dir,
  "dataset_lisi_stratified.tsv.gz"
))
saveRDS(
  umap_plot_data,
  file.path(evaluation_dir, "umap_plot_data.rds"),
  compress = TRUE
)
write_tsv(
  data.frame(
    Metric = "dataset iLISI",
    Space = "method-specific 50-dimensional latent representation",
    Cells = nrow(meta),
    Datasets = length(unique(meta$Dataset)),
    Known_Donors = length(unique(stats::na.omit(known_donor))),
    Perplexity = perplexity,
    Nearest_Neighbor_Error_Bound = nn_eps,
    Numeric_Precision = "double",
    LISI_Source_Commit = unname(lisi_contract[["source_commit"]]),
    LISI_Numeric_Patch = unname(lisi_contract[["numeric_patch"]]),
    LISI_Numeric_Patch_SHA256 =
      unname(lisi_contract[["numeric_patch_sha256"]]),
    LISI_Lower_Bound_Adjustments = sum(
      numeric_audit$Theoretical_Lower_Bound_Adjustments
    ),
    Cell_Order_SHA256 = cell_hash,
    Cell_Order_SHA256_Source = basename(reduction_contract_file),
    Inference_Unit = "source dataset",
    Cell_Level_Use = "descriptive only; no cell-level hypothesis test",
    stringsAsFactors = FALSE
  ),
  file.path(evaluation_dir, "evaluation_contract.tsv")
)

rm(object, meta, lisi_results, umap_plot_data)
gc()
thisutils::log_message("[integration-05] ", 
  paste(
    "Completed full-cell latent-space iLISI and dataset-level uncertainty",
    "for", format(nrow(dataset_scores), big.mark = ","), "datasets"
  )
)
