args <- commandArgs(trailingOnly = TRUE)
overwrite <- length(args) >= 1L &&
  tolower(args[[1L]]) %in% c("true", "t", "1")
metadata_only <- "--metadata-only" %in% args

source("functions/data_paths.R")
source("functions/processed_object.R")

repo_dir <- brainomics_repo_root()
raw_dir <- brainomics_data_path("raw", "GSE294786")
processed_dir <- brainomics_data_path("processed", "GSE294786")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(
  raw_dir,
  "GSE294786_all_counts_transposed.tsv.gz"
)
converted_file <- file.path(
  processed_dir,
  "GSE294786_all_counts_full.h5ad"
)
conversion_audit <- file.path(
  processed_dir,
  "counts_conversion_audit.json"
)
converter <- file.path(
  repo_dir,
  "processing/gse294786_tsv_to_h5ad.py"
)
metadata_file <- file.path(
  raw_dir,
  "GSE294786_all_meta_progenitor.tsv.gz"
)

download_record <- function(file_name) {
  output <- system2(
    "bash",
    c(
      file.path(repo_dir, "download", "GSE294786.sh"),
      "--record",
      file_name
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("cannot resolve GSE294786 download record for ", file_name)
  }
  fields <- strsplit(output[[1L]], "|", fixed = TRUE)[[1L]]
  if (length(fields) != 4L || fields[[2L]] != file_name) {
    stop("invalid GSE294786 download record for ", file_name)
  }
  list(
    bytes = as.numeric(fields[[3L]]),
    sha256 = fields[[4L]]
  )
}

validate_source_file <- function(file, record) {
  if (!file.exists(file)) {
    stop("missing GSE294786 source file: ", file)
  }
  if (as.numeric(file.info(file)$size) != record$bytes) {
    stop("GSE294786 source file size differs from its download record: ", file)
  }
  if (processed_file_sha256(file) != record$sha256) {
    stop("GSE294786 source file SHA-256 differs from its download record: ", file)
  }
}

counts_record <- download_record(basename(counts_file))
metadata_record <- download_record(basename(metadata_file))
validate_source_file(counts_file, counts_record)
validate_source_file(metadata_file, metadata_record)

if (!metadata_only &&
  (overwrite || !all(file.exists(c(converted_file, conversion_audit))))) {
  python <- "python3"
  converter_args <- c(
    converter,
    "--counts", counts_file,
    "--output", converted_file,
    "--audit-output", conversion_audit,
    "--expected-cells", "37994",
    "--expected-features", "23841",
    "--expected-input-sha256", counts_record$sha256
  )
  if (overwrite) {
    converter_args <- c(converter_args, "--overwrite")
  }
  status <- system2(python, shQuote(converter_args))
  if (status != 0L) {
    stop("complete GSE294786 count conversion failed")
  }
}

if (!file.exists(converted_file) || !file.exists(conversion_audit)) {
  stop("complete GSE294786 conversion outputs are missing")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required to validate the GSE294786 conversion audit")
}
audit <- jsonlite::fromJSON(conversion_audit, simplifyVector = TRUE)
if (!isTRUE(audit$all_rows_processed) ||
  !isTRUE(audit$all_features_processed) ||
  audit$input_cells_processed != 37994L ||
  audit$input_features_processed != 23841L ||
  audit$input_sha256 != counts_record$sha256 ||
  audit$output_sha256 != processed_file_sha256(converted_file)) {
  stop("GSE294786 conversion audit differs from the complete source files")
}

PROCESSING_DATASET <- "GSE294786"
source("processing/process_h5ad.R")
