suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")


seed <- as.integer(Sys.getenv(
  "BRAINOMICS_INTEGRATION_SEED",
  unset = "20260730"
))
n_features <- as.integer(Sys.getenv(
  "BRAINOMICS_VARIABLE_FEATURES",
  unset = "3000"
))
expected_reference_cells <- as.numeric(Sys.getenv(
  "BRAINOMICS_EXPECTED_REFERENCE_CELLS",
  unset = "2602031"
))
overwrite <- tolower(Sys.getenv(
  "BRAINOMICS_SCVI_INPUT_OVERWRITE",
  unset = "false"
)) %in% c("true", "t", "1")
res_dir <- brainomics_data_path("integration_25")
input_file <- file.path(res_dir, "objects_filtered.rds")
list_file <- file.path(res_dir, "objects_list_processed.rds")
scvi_input_dir <- file.path(res_dir, "scvi_input")
eligibility_file <- file.path(res_dir, "scvi_input_eligibility.tsv")
checkpoint_dir <- file.path(res_dir, "r_checkpoints")
hvg_checkpoint <- file.path(checkpoint_dir, "objects_hvg.rds")
timing_file <- file.path(res_dir, "integration_stage_timing.tsv")
run_id <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")

record_stage_timing <- function(stage, event, started = NULL) {
  row <- data.frame(
    Run_ID = run_id,
    Stage = stage,
    Event = event,
    Timestamp_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    Elapsed_Seconds = if (is.null(started)) {
      NA_real_
    } else {
      as.numeric(
        difftime(Sys.time(), started, units = "secs")
      )
    },
    stringsAsFactors = FALSE
  )
  write.table(
    row,
    timing_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = !file.exists(timing_file),
    append = file.exists(timing_file)
  )
  invisible(row)
}

if (!file.exists(input_file) || !file.exists(list_file)) {
  stop("Missing filtered reference or processed object-list checkpoint")
}
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

thisutils::log_message("[scvi-input] ", paste("Checking raw counts for all", length(reference_datasets()), "reference datasets"))
stage_started <- Sys.time()
record_stage_timing("scVI_input_audit", "start")
objects_list <- readRDS(list_file)
eligibility <- audit_scvi_input_eligibility(objects_list)
if (nrow(eligibility) != length(reference_datasets()) ||
  any(!eligibility$ScVI_Eligible) ||
  any(!eligibility$Numeric_Count_Eligible) ||
  any(!eligibility$Source_Provenance_Eligible)) {
  stop("Every primary-reference dataset must pass the raw-count scVI gate")
}
write_tsv(eligibility, eligibility_file)
write_tsv(
  eligibility[FALSE, , drop = FALSE],
  file.path(res_dir, "scvi_excluded_datasets.tsv")
)
record_stage_timing("scVI_input_audit", "end", stage_started)

expected_cells <- unlist(lapply(objects_list, colnames), use.names = FALSE)
if (anyDuplicated(expected_cells)) {
  stop("Primary-reference object list contains duplicated cell identifiers")
}
if (length(expected_cells) != expected_reference_cells) {
  stop(
    "Primary-reference cell count differs: expected ",
    format(expected_reference_cells, scientific = FALSE),
    ", observed ",
    length(expected_cells)
  )
}

thisutils::log_message("[scvi-input] ", "Selecting the common feature set before any integration method")
stage_started <- Sys.time()
record_stage_timing("HVG", "start")
objects <- load_processed_object(input_file)
if (!identical(
  sort(unique(as.character(objects$Dataset))),
  sort(reference_datasets())
) || !identical(colnames(objects), expected_cells)) {
  stop("Filtered reference differs from the ordered 22-dataset object list")
}
if (!"Integration_Batch_ID" %in% names(objects[[]]) ||
  anyNA(normalize_missing_metadata(objects$Integration_Batch_ID))) {
  stop("Filtered reference lacks the audited Integration_Batch_ID")
}
count_layers <- Layers(objects, assay = "RNA", search = "^counts")
if (length(count_layers) != length(reference_datasets())) {
  stop("Filtered reference does not retain one count layer per dataset")
}

set.seed(seed)
objects <- NormalizeData(objects, verbose = FALSE)
objects <- FindVariableFeatures(
  objects,
  selection.method = "vst",
  nfeatures = n_features,
  verbose = FALSE
)
selected_features <- VariableFeatures(objects)
if (length(selected_features) != n_features ||
  anyNA(selected_features) ||
  anyDuplicated(selected_features)) {
  stop("Variable-feature selection did not return the fixed feature count")
}
write_tsv(
  data.frame(
    Gene = selected_features,
    HVG_Rank = seq_along(selected_features),
    stringsAsFactors = FALSE
  ),
  file.path(res_dir, "integration_features.tsv")
)
objects@misc$BrainOmics_Integration_Checkpoint <- list(
  Stage = "hvg",
  Seed = seed,
  Dimensions = seq_len(as.integer(Sys.getenv(
    "BRAINOMICS_INTEGRATION_DIMS",
    unset = "50"
  ))),
  Variable_Features = n_features,
  Batch_Model = "dataset",
  Reference_Datasets = reference_datasets(),
  Cells = ncol(objects),
  Features = nrow(objects)
)
save_processed_object(objects, hvg_checkpoint, validate_reload = FALSE)
if (!file.exists(hvg_checkpoint) || file.info(hvg_checkpoint)$size <= 0) {
  stop("HVG checkpoint was not written")
}
record_stage_timing("HVG", "end", stage_started)

thisutils::log_message("[scvi-input] ", "Exporting all 22 datasets for Python scVI")
stage_started <- Sys.time()
record_stage_timing("scVI_export", "start")
input_audit <- export_scvi_input_shards(
  objects_list = objects_list,
  features = selected_features,
  output_dir = scvi_input_dir,
  batch_column = "Integration_Batch_ID",
  eligibility_audit = eligibility,
  expected_cells = expected_cells,
  overwrite = overwrite
)
if (input_audit$Datasets[[1L]] != length(reference_datasets()) ||
  input_audit$Excluded_Datasets[[1L]] != 0L ||
  input_audit$Cells[[1L]] != length(expected_cells)) {
  stop("scVI export does not cover the complete 22-dataset reference")
}
record_stage_timing("scVI_export", "end", stage_started)
write_tsv(
  data.frame(
    Parameter = c(
      "random_seed", "variable_features", "batch_key",
      "reference_datasets", "reference_cells", "excluded_datasets"
    ),
    Value = c(
      seed, n_features, "Integration_Batch_ID",
      length(reference_datasets()), length(expected_cells), 0L
    ),
    stringsAsFactors = FALSE
  ),
  file.path(res_dir, "scvi_input_parameters.tsv")
)

rm(objects, objects_list)
gc()
thisutils::log_message("[scvi-input] ", paste(
  "Saved complete scVI input:",
  format(length(expected_cells), big.mark = ","),
  "cells from",
  length(reference_datasets()),
  "raw-count datasets"
))
