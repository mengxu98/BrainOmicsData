#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) {
    return(default)
  }
  args[[index + 1L]]
}

package_dir <- normalizePath(
  value_after("--package-dir", "../../data/BrainOmicsData/ScienceDB"),
  mustWork = TRUE
)
output_dir <- value_after(
  "--output-dir", "results/integration_25/publication_tables"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- file.path(package_dir, "provenance", "file_manifest.tsv")
md5_file <- file.path(package_dir, "md5sum.txt")
required <- c(manifest_file, md5_file)
if (any(!file.exists(required))) {
  stop("Final release manifest or checksum chains are missing")
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required to validate the final package")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table is required to validate the public metadata")
}

manifest <- utils::read.delim(
  manifest_file,
  sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)
required_columns <- c(
  "path", "file_name", "size_bytes", "rows", "columns",
  "nonzero_values", "schema_version", "license_scope", "description",
  "md5", "sha256"
)
if (!identical(names(manifest), required_columns) ||
  anyDuplicated(manifest$path) || any(grepl("^/", manifest$path))) {
  stop("Final package manifest schema or relative-path contract changed")
}
payload_files <- file.path(package_dir, manifest$path)
if (any(!file.exists(payload_files)) || any(file.info(payload_files)$isdir)) {
  stop("Final package manifest contains a missing payload file")
}
observed_size <- as.numeric(file.info(payload_files)$size)
observed_md5 <- unname(tools::md5sum(payload_files))
observed_sha256 <- vapply(
  payload_files,
  digest::digest,
  character(1L),
  algo = "sha256", file = TRUE, serialize = FALSE
)
if (!identical(as.numeric(manifest$size_bytes), observed_size) ||
  !identical(manifest$md5, observed_md5) ||
  !identical(manifest$sha256, unname(observed_sha256))) {
  stop("Final package payload differs from its post-anonymization manifest")
}

metadata_file <- file.path(package_dir, "metadata", "metadata.tsv.gz")
public_metadata <- data.table::fread(
  metadata_file,
  sep = "\t", data.table = FALSE,
  na.strings = c("", "NA", "Not reported")
)
expected_metadata_columns <- c(
  "Cells", "Dataset", "Original_Sample_ID", "Donor_ID", "Sample_ID",
  "Library_ID", "Technology", "Modality", "Age", "AgeIntervalID", "Sex",
  "BrainRegion", "Source_CellType", "Cluster", "CellType"
)
if (!identical(names(public_metadata), expected_metadata_columns) ||
  anyDuplicated(public_metadata$Cells) ||
  any(!grepl("^Cell[0-9]{7}$", public_metadata$Cells))) {
  stop("Public metadata schema, cell identity or cell uniqueness changed")
}
missing_public_value <- function(values) {
  values <- trimws(as.character(values))
  is.na(values) | !nzchar(values) | values %in% c("NA", "Not reported")
}
technology_missing <- sum(missing_public_value(public_metadata$Technology))
modality_missing <- sum(missing_public_value(public_metadata$Modality))
original_sample_id_missing <- sum(missing_public_value(
  public_metadata$Original_Sample_ID
))
if (original_sample_id_missing != 0L ||
  technology_missing != 0L || modality_missing != 0L) {
  stop(
    "Public metadata coverage is incomplete: ",
    "Original_Sample_ID=", original_sample_id_missing,
    "; Technology=", technology_missing, "; Modality=", modality_missing
  )
}
metadata_manifest_rows <- manifest$rows[
  manifest$path == "metadata/metadata.tsv.gz"
]
if (length(metadata_manifest_rows) != 1L ||
  as.numeric(metadata_manifest_rows) != nrow(public_metadata)) {
  stop("Public metadata row count differs from the release manifest")
}

component <- sub("/.*$", "", manifest$path)
component[!grepl("/", manifest$path, fixed = TRUE)] <- "package root"
component_levels <- c(
  "expression", "metadata", "embeddings", "objects", "validation",
  "provenance", "scripts", "package root"
)
unknown <- setdiff(unique(component), component_levels)
if (length(unknown) > 0L) {
  stop("Unregistered package component: ", paste(unknown, collapse = ", "))
}
component_description <- c(
  expression = "Complete dataset-resolved raw-count shards and ordered shard manifest",
  metadata = "Cell, feature, dataset and schema metadata",
  embeddings = "Cell-ordered integrated and unintegrated coordinate tables",
  objects = "Lightweight Seurat marker-plotting and metadata object",
  validation = "Integration, annotation, lineage and held-out validation outputs",
  provenance = "Source access, crosswalk, inclusion and audit records",
  scripts = "R and Python reader examples",
  `package root` = "README and package-level documentation"
)
summary_rows <- lapply(component_levels, function(current) {
  selected <- component == current
  if (!any(selected)) {
    return(NULL)
  }
  data.frame(
    Component = current,
    Payload_Files = sum(selected),
    Size_Bytes = sum(observed_size[selected]),
    Size_GiB = sum(observed_size[selected]) / 1024^3,
    Description = unname(component_description[[current]]),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
if (sum(summary$Payload_Files) != nrow(manifest) ||
  sum(summary$Size_Bytes) != sum(observed_size)) {
  stop("Package component summary does not cover every payload exactly once")
}

chain_files <- c(manifest_file, md5_file)
contract <- data.frame(
  Field = c(
    "Payload_Files", "Manifest_and_Checksum_Files", "Total_Files",
    "Payload_Size_Bytes", "Total_Package_Size_Bytes",
    "Manifest_Schema_Version", "Public_Metadata_Rows",
    "Public_Metadata_Columns", "Original_Sample_ID_Missing_Cells",
    "Technology_Missing_Cells", "Modality_Missing_Cells",
    "Validation_Status"
  ),
  Value = c(
    nrow(manifest), length(chain_files), nrow(manifest) + length(chain_files),
    sum(observed_size),
    sum(observed_size) + sum(as.numeric(file.info(chain_files)$size)),
    "2.0.0", nrow(public_metadata), ncol(public_metadata),
    original_sample_id_missing, technology_missing, modality_missing,
    "complete post-anonymization payload validation"
  ),
  stringsAsFactors = FALSE
)

summary_file <- file.path(output_dir, "sciencedb_package_summary.tsv")
contract_file <- file.path(output_dir, "sciencedb_package_contract.tsv")
utils::write.table(
  summary, summary_file,
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)
utils::write.table(
  contract, contract_file,
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)
cat(
  "[sciencedb-package-summary] PASS:", nrow(manifest), "payload files;",
  nrow(manifest) + length(chain_files), "total files;",
  sprintf("%.3f GiB", sum(as.numeric(file.info(c(payload_files, chain_files))$size)) / 1024^3),
  "\n"
)
