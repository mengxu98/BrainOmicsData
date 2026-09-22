#!/usr/bin/env Rscript

source("functions/data_paths.R")
source("functions/metadata_schema.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/original_source_readers.R")
source("functions/integration.R")

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (exists("PROCESSING_DATASET", inherits = FALSE)) {
  get("PROCESSING_DATASET", inherits = FALSE)
} else if (length(args) >= 1L) {
  args[[1L]]
} else {
  stop("Usage: Rscript processing/process_original_source.R <dataset> [overwrite]")
}
overwrite_arg <- if (exists("PROCESSING_DATASET", inherits = FALSE)) {
  if (length(args) >= 1L) args[[1L]] else "false"
} else {
  if (length(args) >= 2L) args[[2L]] else "false"
}
overwrite <- tolower(overwrite_arg) %in% c("true", "t", "1")
supported <- c(
  "AllenM1", "EGAS00001006537",
)
if (!dataset %in% supported) {
  stop("Unsupported original-source dataset: ", dataset)
}

require_original_source_packages()
data_root <- brainomics_data_root()
raw_dir <- brainomics_data_path("raw", dataset)
processed_dir <- brainomics_data_path("processed", dataset)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(processed_dir, paste0(dataset, "_processed.rds"))
audit_file <- file.path(processed_dir, "source_reconstruction_audit.tsv")
if (file.exists(output_file) && !overwrite) {
  message("Retaining existing complete original-source object: ", output_file)
  quit(save = "no", status = 0L)
}

source("processing/original_source_common.R")
source("processing/original_source_geo10x.R")
source("processing/AllenM1.R")
source("processing/EGAS00001006537.R")

message("Reconstructing complete original source for ", dataset)
object <- switch(dataset,
  AllenM1 = process_allen(),
  EGAS00001006537 = process_egas(),
  GSE202210 = process_gse202210()
)
if (!identical(colnames(object), rownames(object[[]]))) {
  object <- object[, rownames(object[[]])]
}
validate_processed_object(object)
save_processed_object(object, output_file)
write_source_reconstruction_audit(object, audit_file, dataset)
message("Saved complete original-source object: ", output_file)
