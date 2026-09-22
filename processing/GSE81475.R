source("functions/prepare_env.R")
source("functions/data_paths.R")
source("functions/dataset_metadata.R")

data_dir <- brainomics_data_path("raw", "GSE81475")
res_dir <- check_dir(brainomics_data_path("processed", "GSE81475"))
counts_file <- file.path(data_dir, "counts.txt")
metadata_file <- file.path(data_dir, "metadata.csv")
series_file <- file.path(data_dir, "GSE81475_series_matrix.txt.gz")

for (file in c(counts_file, metadata_file, series_file)) {
  if (!file.exists(file)) {
    stop("GSE81475 required source file is missing: ", file)
  }
}

thisutils::log_message("Loading the complete GSE81475 count table and GEO metadata...")
counts_source <- read.delim(
  counts_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
if (ncol(counts_source) != 1608L) {
  stop("GSE81475 complete GEO count table must contain 1,608 cells")
}
metadata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (nrow(metadata) != 476L) {
  stop("GSE81475 retained fetal DFC metadata must contain 476 cells")
}
geo <- read_geo_series_matrix_sample_metadata(series_file)
required_geo <- c(
  "Sample_geo_accession", "Sample_title", "Sample_instrument_model",
  "characteristic_donor", "characteristic_age",
  "characteristic_c1_chip"
)
missing_geo <- setdiff(required_geo, names(geo))
if (length(missing_geo) > 0L || nrow(geo) != 1608L) {
  stop(
    "GSE81475 GEO series metadata is incomplete: ",
    paste(missing_geo, collapse = ", ")
  )
}

normalize_cell_id <- function(value) {
  gsub("[^A-Za-z0-9]", "", as.character(value))
}
geo_index <- match(
  normalize_cell_id(rownames(metadata)),
  normalize_cell_id(geo$Sample_title)
)
if (anyNA(geo_index) || anyDuplicated(geo_index)) {
  stop("GSE81475 retained cells do not map one-to-one to GEO records")
}
geo_retained <- geo[geo_index, , drop = FALSE]
count_index <- match(
  normalize_cell_id(geo_retained$Sample_title),
  normalize_cell_id(colnames(counts_source))
)
if (anyNA(count_index) || anyDuplicated(count_index)) {
  stop("GSE81475 retained GEO records do not map to the count matrix")
}

donor <- as.character(geo_retained$characteristic_donor)
age <- as.character(geo_retained$characteristic_age)
capture_device <- as.character(geo_retained$characteristic_c1_chip)
sequencing_lane <- sub(
  ".*_(L[0-9]{3})_R[12]_.*$",
  "\\1",
  geo_retained$Sample_title
)
if (
  length(unique(donor)) != 3L ||
    length(unique(capture_device)) != 6L ||
    length(unique(geo_retained$Sample_geo_accession)) != 476L ||
    any(!grepl("^L[0-9]{3}$", sequencing_lane))
) {
  stop("GSE81475 audited 3-donor/6-chip/476-GSM design differs")
}
chip_donors <- vapply(
  split(donor, capture_device),
  function(value) length(unique(value)),
  integer(1)
)
if (any(chip_donors != 1L)) {
  stop("GSE81475 each Fluidigm C1 chip must remain within one donor")
}

thisutils::log_message("Aggregating duplicated gene symbols by count summation...")
retained_counts <- as.matrix(
  counts_source[, count_index, drop = FALSE]
)
gene_names <- sub("^[^|]+[|]", "", rownames(retained_counts))
if (anyNA(gene_names) || any(gene_names == "")) {
  stop("GSE81475 contains features without a gene symbol")
}
counts <- rowsum(
  retained_counts,
  group = gene_names,
  reorder = FALSE
)
counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
colnames(counts) <- geo_retained$Sample_title
rm(counts_source, retained_counts)
gc()

metadata <- metadata[geo_index * 0L + seq_len(nrow(metadata)), , drop = FALSE]
rownames(metadata) <- geo_retained$Sample_title
metadata$Cells <- geo_retained$Sample_title
metadata$Dataset <- "GSE81475"
metadata$Technology <- "Fluidigm C1"
metadata$Sequence <- "scRNA-seq"
metadata$Assay_Type <- "single-cell RNA sequencing"
metadata$Sequencing_Platform <- geo_retained$Sample_instrument_model
metadata$Library_Chemistry <- paste(
  "Fluidigm C1 capture; SMARTer Ultra Low RNA cDNA;",
  "Illumina Nextera XT library"
)
metadata$Original_Donor_ID <- donor
metadata$Original_Specimen_ID <- paste(donor, "DFC", sep = "::")
metadata$Original_Library_ID <- geo_retained$Sample_title
metadata$Original_Source_Sample_ID <- paste(donor, "DFC", sep = "::")
metadata$Original_Source_Record_ID <- geo_retained$Sample_geo_accession
metadata$Original_Technical_Batch_ID <- capture_device
metadata$Source_GEO_Record_ID <- geo_retained$Sample_geo_accession
metadata$Source_Single_Cell_Library_ID <- geo_retained$Sample_title
metadata$Source_Capture_Device_ID <- capture_device
metadata$Source_Sequencing_Lane <- sequencing_lane
metadata$Source_Full_Series_Cell_Count <- 1608L
metadata$Source_Retained_Cell_Count <- 476L
metadata$Source_Selection_Rule <- paste(
  "retain the existing non-infected fetal DFC reference subset;",
  "exclude neural epithelial stem-cell culture records and other",
  "source cells not present in the original 476-cell reference metadata"
)
metadata$Sample <- metadata$Original_Specimen_ID
metadata$Sample_ID <- metadata$Original_Source_Record_ID
metadata$Brain_Region <- "Dorsolateral prefrontal cortex"
metadata$Region <- "DFC"
metadata$CellType_raw <- as.character(metadata$CellType)
metadata <- standardize_source_age_metadata(
  metadata,
  source_age = age,
  source_basis = "postconceptional",
  source_unit = "post-conception weeks",
  source_reference = paste(
    "GEO GSE81475 characteristics and Li et al. 2016,",
    "DOI 10.1016/j.celrep.2016.08.038"
  )
)
metadata$Sex_Source_Raw <- as.character(metadata$Sex)
metadata$Sex_Assignment_Method <- "source study metadata"
metadata$Sex <- standardize_source_sex_value(metadata$Sex_Source_Raw)

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Assay_Type",
  "Sequencing_Platform", "Library_Chemistry", "Sample", "Sample_ID",
  "Original_Donor_ID", "Original_Specimen_ID", "Original_Library_ID",
  "Original_Source_Sample_ID", "Original_Source_Record_ID",
  "Original_Technical_Batch_ID", "Source_GEO_Record_ID",
  "Source_Single_Cell_Library_ID", "Source_Capture_Device_ID",
  "Source_Sequencing_Lane", "CellType_raw", "Brain_Region", "Region",
  "Age_Source_Raw", "Age_Source_Unit", "Age_Source_Basis",
  "Age_Source_Reference", "Age_Harmonization_Input",
  "Age_Conversion_Formula", "Age_Conversion_Confidence",
  "Age_Conversion_Applied", "Age", "Sex_Source_Raw",
  "Sex_Assignment_Method", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)
counts <- counts[, metadata$Cells, drop = FALSE]

thisutils::log_message("Creating the complete retained GSE81475 Seurat object...")
object <- CreateSeuratObject(counts = counts, meta.data = metadata)
saveRDS(object, file.path(res_dir, "GSE81475_processed.rds"))
