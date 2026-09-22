suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(future)
})

source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")


future_configuration <- configure_complete_integration_future()

seed <- as.integer(Sys.getenv(
  "BRAINOMICS_INTEGRATION_SEED",
  unset = "20260730"
))
dims <- seq_len(as.integer(Sys.getenv(
  "BRAINOMICS_INTEGRATION_DIMS",
  unset = "50"
)))
remove_intermediates <- tolower(Sys.getenv(
  "BRAINOMICS_REMOVE_INTEGRATION_INTERMEDIATES",
  unset = "false"
)) %in% c("true", "t", "1")
res_dir <- brainomics_data_path("integration_25")
input_file <- file.path(res_dir, "objects_integrated_r.rds")
output_file <- file.path(res_dir, "objects_integrated.rds")
scvi_output_dir <- file.path(res_dir, "scvi_output")
eligibility_file <- file.path(res_dir, "scvi_input_eligibility.tsv")
if (!file.exists(input_file)) {
  stop("Missing R-integrated reference object: ", input_file)
}
if (!file.exists(eligibility_file)) {
  stop("Missing scVI input eligibility audit: ", eligibility_file)
}

thisutils::log_message("[integration-04] ", paste(
  "Loading the R-integrated",
  length(reference_datasets()),
  "-dataset reference"
))
objects <- load_processed_object(input_file)
if (!identical(
  sort(unique(as.character(objects$Dataset))),
  sort(reference_datasets())
)) {
  stop("R-integrated reference does not contain the expected datasets")
}

scvi_eligibility <- read.delim(
  eligibility_file,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!identical(as.character(scvi_eligibility$Dataset), reference_datasets()) ||
  anyNA(scvi_eligibility$ScVI_Eligible) ||
  any(!scvi_eligibility$ScVI_Eligible)) {
  stop("scVI eligibility audit differs from the reference contract")
}
scvi_datasets <- reference_datasets()
scvi_cells <- colnames(objects)
if (length(scvi_cells) == 0L || anyDuplicated(scvi_cells)) {
  stop("scVI reference cohort is empty or duplicated")
}

thisutils::log_message("[integration-04] ", "Importing the source-verified Python scVI latent representation")
scvi_latent <- read_scvi_latent(
  output_dir = scvi_output_dir,
  expected_cells = scvi_cells,
  expected_dims = length(dims)
)
if (!identical(rownames(scvi_latent), scvi_cells)) {
  stop("scVI latent cell order differs after reference reordering")
}
objects[["integrated.scvi"]] <- CreateDimReducObject(
  embeddings = scvi_latent,
  key = "scVI_",
  assay = DefaultAssay(objects)
)
rm(scvi_latent)
gc()

set.seed(seed)
objects <- RunUMAP(
  objects,
  reduction = "integrated.scvi",
  dims = dims,
  reduction.name = "umap.scvi",
  seed.use = seed,
  verbose = FALSE
)

method_identity <- integration_method_identity()
missing_reductions <- setdiff(
  method_identity$Reduction,
  names(objects@reductions)
)
if (length(missing_reductions) > 0L) {
  stop(
    "Final integration reductions are missing: ",
    paste(missing_reductions, collapse = ", ")
  )
}
reduction_audit <- do.call(
  rbind,
  lapply(seq_len(nrow(method_identity)), function(index) {
    reduction <- method_identity$Reduction[[index]]
    embedding <- Embeddings(objects, reduction = reduction)
    expected_reduction_cells <- if (
      identical(method_identity$Display_Name[[index]], "scVI")
    ) {
      scvi_cells
    } else {
      colnames(objects)
    }
    if (!identical(rownames(embedding), expected_reduction_cells) ||
      any(!is.finite(embedding))) {
      stop(reduction, " embedding failed cell-order or finite-value checks")
    }
    data.frame(
      Output_Name = method_identity$Output_Name[[index]],
      Display_Name = method_identity$Display_Name[[index]],
      Reduction = reduction,
      Space = method_identity$Space[[index]],
      Cells = nrow(embedding),
      Dimensions = ncol(embedding),
      Cell_Order_SHA256 = sha256_text_vector(rownames(embedding)),
      Finite = TRUE,
      stringsAsFactors = FALSE
    )
  })
)
comparison_cohort <- do.call(rbind, lapply(
  reference_datasets(),
  function(dataset) {
    dataset_cells <- sum(as.character(objects$Dataset) == dataset)
    data.frame(
      Dataset = dataset,
      Reference_Cells = dataset_cells,
      Four_Method_Comparison_Cells = dataset_cells,
      ScVI_Eligible = TRUE,
      Comparison_Role = "Raw-scVI-Harmony-RPCA common-input comparison",
      stringsAsFactors = FALSE
    )
  }
))
if (sum(comparison_cohort$Four_Method_Comparison_Cells) !=
  length(scvi_cells)) {
  stop("four-method comparison cohort cell count differs")
}
write_tsv(
  method_identity,
  file.path(res_dir, "integration_method_identity.tsv")
)
write_tsv(
  comparison_cohort,
  file.path(res_dir, "integration_comparison_cohort.tsv")
)
write_tsv(
  reduction_audit,
  file.path(res_dir, "integration_reduction_audit.tsv")
)

expected_cells <- ncol(objects)
expected_features <- nrow(objects)
save_processed_object(objects, output_file, validate_reload = FALSE)
if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
  stop("Final integrated reference object was not written")
}
rm(objects)
gc()
validated <- load_processed_object(output_file)
validate_processed_object(
  validated,
  expected_features = expected_features,
  expected_cells = expected_cells
)
if (!identical(
  sort(unique(as.character(validated$Dataset))),
  sort(reference_datasets())
)) {
  stop("Reloaded final reference does not contain the expected datasets")
}
missing_reductions <- setdiff(
  method_identity$Reduction,
  names(validated@reductions)
)
if (length(missing_reductions) > 0L) {
  stop(
    "Reloaded final integration reductions are missing: ",
    paste(missing_reductions, collapse = ", ")
  )
}
if (!identical(
  rownames(Embeddings(validated, reduction = "integrated.scvi")),
  scvi_cells
)) {
  stop("Reloaded scVI reduction differs from the common-cell cohort")
}
rm(validated)
gc()
if (remove_intermediates) {
  unlink(input_file)
  if (file.exists(input_file)) {
    stop("Could not remove the superseded R integration checkpoint")
  }
}
thisutils::log_message("[integration-04] ", 
  paste(
    "Saved the common-input",
    length(reference_datasets()),
    "-dataset four-method reference:",
    format(expected_cells, big.mark = ","),
    "cells"
  )
)
