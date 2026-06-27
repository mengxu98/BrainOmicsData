#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

out_dir <- normalizePath(
  value_after("--out-dir", "../../data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource"),
  mustWork = TRUE
)
readme_dir <- file.path(out_dir, "00_README")
dir.create(readme_dir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- file.path(readme_dir, "file_manifest.tsv")
dictionary_file <- file.path(readme_dir, "data_dictionary.tsv")

write_tsv <- function(x, file) {
  utils::write.table(
    x, file = file, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "missing"
  )
}

relative_path <- function(files) {
  sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", out_dir), "/?"), "", files)
}

table_connection <- function(file) {
  if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
}

table_shape <- function(file) {
  if (!grepl("\\.(tsv|tsv\\.gz|csv|csv\\.gz)$", file)) {
    return(c(rows = "not applicable", columns = "not applicable"))
  }
  con <- table_connection(file)
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1, warn = FALSE)
  if (length(header) == 0) return(c(rows = "0", columns = "0"))
  sep <- if (grepl("\\.csv(\\.gz)?$", file)) "," else "\t"
  columns <- length(strsplit(header, sep, fixed = TRUE)[[1]])
  rows <- 0L
  repeat {
    chunk <- readLines(con, n = 100000, warn = FALSE)
    if (length(chunk) == 0) break
    rows <- rows + length(chunk)
  }
  c(rows = as.character(rows), columns = as.character(columns))
}

read_header_example <- function(file) {
  if (!grepl("\\.(tsv|tsv\\.gz|csv|csv\\.gz)$", file)) return(NULL)
  con <- table_connection(file)
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, n = 2, warn = FALSE)
  if (length(lines) == 0) return(NULL)
  sep <- if (grepl("\\.csv(\\.gz)?$", file)) "," else "\t"
  header <- strsplit(lines[[1]], sep, fixed = TRUE)[[1]]
  example <- if (length(lines) >= 2) strsplit(lines[[2]], sep, fixed = TRUE)[[1]] else rep("missing", length(header))
  length(example) <- length(header)
  data.frame(column_name = header, example_value = example, stringsAsFactors = FALSE)
}

folder_section <- function(folder) {
  if (startsWith(folder, "00_README")) return("README and package documentation")
  if (startsWith(folder, "01_source_datasets")) return("Data Records - source datasets")
  if (startsWith(folder, "02_sample_metadata")) return("Data Records - sample metadata")
  if (startsWith(folder, "03_cell_metadata")) return("Data Records - cell metadata")
  if (startsWith(folder, "04_processed_objects")) return("Data Records - processed objects")
  if (startsWith(folder, "05_embeddings_and_clusters")) return("Technical Validation - embeddings and clusters")
  if (startsWith(folder, "06_quality_control")) return("Technical Validation - quality control")
  if (startsWith(folder, "07_celltype_annotation")) return("Technical Validation - cell-type annotation")
  if (startsWith(folder, "08_reusable_workflow")) return("Methods - reusable workflow")
  if (startsWith(folder, "09_value_added_annotations")) return("Usage Notes - optional annotations")
  if (startsWith(folder, "10_figures_and_tables")) return("Figure and table source data")
  "Data Records"
}

access_status <- function(path) {
  if (grepl("controlled_access|access_level", path)) return("metadata_only_for_restricted_sources")
  if (grepl("04_processed_objects", path)) return("processed_or_derived_data")
  if (grepl("09_value_added_annotations", path)) return("optional_derived_annotation")
  "open_metadata_or_derived_output"
}

describe_file <- function(path) {
  name <- basename(path)
  if (name == "README.md") return("Package-level README and reuse notes.")
  if (name == "file_manifest.tsv") return("Machine-readable inventory of all files in the ScienceDB package.")
  if (name == "data_dictionary.tsv") return("Field-level data dictionary in TSV format.")
  if (name == "data_dictionary.xlsx") return("Field-level data dictionary in Excel format.")
  if (name == "source_dataset_summary.tsv") return("Summary of source datasets, accessions, repositories, access levels and retained cells.")
  if (name == "source_accessions.tsv") return("Original repository accessions and access routes for source data.")
  if (name == "sample_metadata_harmonized.tsv") return("Harmonized sample-level metadata with missing fields explicitly marked.")
  if (name == "cell_metadata_harmonized.tsv.gz") return("Harmonized cell-level metadata for retained cells or nuclei.")
  if (name == "age_interval_mapping.tsv") return("Mapping between S1-S15 age interval IDs, interval labels and age ranges.")
  if (name == "brain_region_mapping.tsv") return("Mapping of harmonized brain-region labels and optional region hierarchy.")
  if (name == "celltype_cluster_mapping.tsv") return("Mapping from integrated Seurat cluster IDs to major cell-type annotations.")
  if (name == "umap_coordinates.tsv.gz") return("UMAP coordinates for retained cells or nuclei.")
  if (name == "pca_coordinates.tsv.gz") return("PCA or integrated RPCA coordinates for retained cells or nuclei.")
  if (name == "qc_metrics_by_sample.tsv") return("Sample-level quality-control and coverage metrics.")
  if (name == "qc_metrics_by_dataset.tsv") return("Dataset-level quality-control and coverage metrics.")
  if (name == "lisi_before_after_integration.tsv") return("Dataset-label LISI summaries before and after integration.")
  if (name == "marker_gene_list.tsv") return("Marker genes used for manual major cell-type annotation.")
  if (name == "marker_expression_summary.tsv") return("Marker expression summaries used to validate cell-type labels.")
  if (name == "controlled_access_notes.txt") return("Notes for controlled-access raw data and redistribution limits.")
  paste("ScienceDB package file:", name)
}

format_for <- function(path) {
  if (grepl("\\.tsv\\.gz$", path)) return("TSV.GZ")
  if (grepl("\\.tsv$", path)) return("TSV")
  if (grepl("\\.csv\\.gz$", path)) return("CSV.GZ")
  if (grepl("\\.csv$", path)) return("CSV")
  if (grepl("\\.xlsx$", path)) return("XLSX")
  if (grepl("\\.md$", path)) return("Markdown")
  if (grepl("\\.txt$", path)) return("Text")
  if (grepl("\\.rds$", path)) return("RDS")
  if (grepl("\\.h5ad$", path)) return("H5AD")
  if (grepl("\\.h5seurat$", path)) return("H5Seurat")
  tools::file_ext(path)
}

column_description <- function(column_name) {
  descriptions <- c(
    cell_id = "Unique cell or nucleus identifier in the harmonized resource.",
    source_cell_id = "Original cell barcode or source identifier from the input dataset.",
    sample_id = "Harmonized sample identifier used across sample and cell metadata.",
    donor_id = "Anonymized donor or donor-like grouping identifier where available.",
    source_dataset = "Dataset label used by the BrainOmicsData/HARNexus workflow.",
    source_accession = "Original accession or stable identifier for the source data.",
    source_repository = "Repository or database where the source data are hosted.",
    publication_reference = "Original publication DOI or reference string.",
    publication_doi = "Digital Object Identifier for the source publication, where resolved.",
    access_level = "Access category: public, controlled, restricted or mixed.",
    brain_region_original = "Original brain-region label as reported by the source, if retained.",
    brain_region_harmonized = "Standardized brain-region label used in this resource.",
    region_level_1 = "Optional broad brain-region hierarchy.",
    region_level_2 = "Optional finer brain-region hierarchy.",
    reported_age = "Age value as reported by the source dataset.",
    reported_age_unit = "Reported age unit after harmonization, such as PCW or years.",
    age_interval_id = "Unified S1-S15 age interval identifier used for integration and querying.",
    age_interval = "Human-readable label for the unified age interval.",
    age_range = "Age range covered by the unified age interval.",
    sex = "Reported donor sex, or missing/not reported.",
    sequencing_modality = "Sequencing modality or source assay description.",
    single_cell_or_single_nucleus = "Whether the measurement is single-cell or single-nucleus.",
    sequencing_platform = "Sequencing platform or technology label.",
    library_protocol = "Library preparation protocol when reported.",
    cluster_id = "Integrated Seurat cluster identifier.",
    cell_type = "Major marker-based cell-type annotation.",
    n_counts = "Total molecule/count value per cell where available.",
    n_genes = "Detected gene/feature count per cell where available.",
    percent_mito = "Percentage of mitochondrial counts per cell where available.",
    percent_ribo = "Percentage of ribosomal counts per cell where available.",
    quality_control_status = "Whether the record passed the exported QC/metadata filter.",
    cell_count = "Number of retained cells or nuclei.",
    marker_gene = "Marker gene used for cell-type annotation or validation.",
    method = "Integration or embedding method.",
    median = "Median metric value.",
    q25 = "First quartile.",
    q75 = "Third quartile.",
    mean = "Mean metric value.",
    max = "Maximum metric value.",
    md5 = "MD5 checksum."
  )
  if (column_name %in% names(descriptions)) descriptions[[column_name]] else paste("Column exported from", column_name, "in the ScienceDB workflow.")
}

data_type <- function(example) {
  if (is.na(example) || example == "" || example == "missing") return("string")
  if (grepl("^-?[0-9]+$", example)) return("integer")
  if (grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", example)) return("numeric")
  "string"
}

files <- list.files(out_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
files <- files[file.info(files)$isdir == FALSE]
rel <- relative_path(files)

shapes <- lapply(files, table_shape)
manifest <- data.frame(
  file_name = basename(rel),
  folder = dirname(rel),
  description = vapply(rel, describe_file, character(1)),
  format = vapply(rel, format_for, character(1)),
  rows = vapply(shapes, function(x) x[["rows"]], character(1)),
  columns = vapply(shapes, function(x) x[["columns"]], character(1)),
  compressed = ifelse(grepl("\\.gz$", rel), "yes", "no"),
  md5 = unname(tools::md5sum(files)),
  access_status = vapply(rel, access_status, character(1)),
  related_section = vapply(dirname(rel), folder_section, character(1)),
  path = rel,
  size_bytes = as.numeric(file.info(files)$size),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$folder, manifest$file_name), ]
write_tsv(manifest, manifest_file)

dict <- do.call(rbind, lapply(seq_along(files), function(i) {
  example <- read_header_example(files[[i]])
  if (is.null(example)) return(NULL)
  rel_i <- rel[[i]]
  data.frame(
    file_name = rel_i,
    column_name = example$column_name,
    description = vapply(example$column_name, column_description, character(1)),
    data_type = vapply(example$example_value, data_type, character(1)),
    allowed_values = "not enumerated",
    missing_value_definition = "missing or not reported indicates source field unavailable; not applicable indicates field does not apply.",
    example_value = ifelse(is.na(example$example_value) | example$example_value == "", "missing", example$example_value),
    file = rel_i,
    field = example$column_name,
    definition = vapply(example$column_name, column_description, character(1)),
    stringsAsFactors = FALSE
  )
}))
if (!is.null(dict) && nrow(dict) > 0) {
  dict <- dict[order(dict$file_name, dict$column_name), ]
  write_tsv(dict, dictionary_file)
}

message("Manifest written: ", manifest_file)
message("Data dictionary TSV written: ", dictionary_file)
