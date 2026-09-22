source("functions/prepare_env.R")
source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")

raw_dir <- brainomics_data_path("raw", "Li_et_al_2018", "rawData")
res_dir <- check_dir(brainomics_data_path("processed", "Li_et_al_2018"))
fetal_file <- file.path(
  raw_dir,
  "Prenatal_scRNA_seq",
  "Sestan.fetalHuman.Psychencode.Rdata"
)
adult_file <- file.path(
  raw_dir,
  "Adult_snRNA_seq",
  "Sestan.adultHumanNuclei.Psychencode.Rdata"
)
gse_series_file <- brainomics_data_path(
  "raw", "GSE81475", "GSE81475_series_matrix.txt.gz"
)
for (file in c(fetal_file, adult_file, gse_series_file)) {
  if (!file.exists(file)) {
    stop("Li et al. 2018 required source file is missing: ", file)
  }
}

load_source_object <- function(file, required) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(file, envir = environment)
  missing <- setdiff(required, loaded)
  if (length(missing) > 0L) {
    stop("Li et al. source RData is missing: ", paste(missing, collapse = ", "))
  }
  mget(required, envir = environment, inherits = FALSE)
}

aggregate_gene_symbols <- function(counts) {
  counts <- as.matrix(counts)
  if (anyNA(counts) || any(counts < 0) || any(counts != floor(counts))) {
    stop("Li et al. source count matrix contains invalid counts")
  }
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  symbols <- sub("^[^|]+[|]", "", rownames(counts))
  if (anyNA(symbols) || any(symbols == "")) {
    stop("Li et al. source matrix contains features without gene symbols")
  }
  levels <- unique(symbols)
  aggregator <- Matrix::sparseMatrix(
    i = match(symbols, levels),
    j = seq_along(symbols),
    x = 1,
    dims = c(length(levels), length(symbols)),
    dimnames = list(levels, rownames(counts))
  )
  result <- methods::as(aggregator %*% counts, "dgCMatrix")
  colnames(result) <- colnames(counts)
  result
}

align_sparse_features <- function(counts, features) {
  missing <- setdiff(features, rownames(counts))
  if (length(missing) > 0L) {
    zeros <- Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      dims = c(length(missing), ncol(counts)),
      dimnames = list(missing, colnames(counts))
    )
    counts <- rbind(counts, zeros)
  }
  methods::as(counts[features, , drop = FALSE], "dgCMatrix")
}

thisutils::log_message("Loading the original Li et al. 2018 fetal and adult RData...")
fetal <- load_source_object(fetal_file, c("count2", "meta2"))
adult <- load_source_object(adult_file, c("umi2", "umi.raw", "meta2"))
if (
  !identical(rownames(adult$umi.raw), rownames(adult$umi2)) ||
    !identical(colnames(adult$umi2), rownames(adult$meta2)) ||
    !all(colnames(adult$umi2) %in% colnames(adult$umi.raw)) ||
    ncol(adult$umi.raw) != 17335L ||
    ncol(adult$umi2) != 17093L
) {
  stop(
    "Li adult raw/retained matrices must preserve the reviewed ",
    "17,335-to-17,093 nucleus relationship"
  )
}
fetal_counts <- aggregate_gene_symbols(fetal$count2)
adult_counts <- aggregate_gene_symbols(
  adult$umi.raw[, colnames(adult$umi2), drop = FALSE]
)
fetal_meta <- fetal$meta2
adult_meta <- adult$meta2
rm(fetal, adult)
gc()

if (
  ncol(fetal_counts) != 762L || nrow(fetal_meta) != 762L ||
    ncol(adult_counts) != 17093L || nrow(adult_meta) != 17093L ||
    !identical(colnames(fetal_counts), rownames(fetal_meta)) ||
    !identical(colnames(adult_counts), rownames(adult_meta))
) {
  stop("Li et al. audited 762-fetal/17,093-adult source design differs")
}

geo <- read_geo_series_matrix_sample_metadata(gse_series_file)
normalize_cell_id <- function(value) {
  gsub("[^A-Za-z0-9]", "", as.character(value))
}
geo_index <- match(
  normalize_cell_id(colnames(fetal_counts)),
  normalize_cell_id(geo$Sample_title)
)
if (anyNA(geo_index) || anyDuplicated(geo_index)) {
  stop("all 762 Li fetal cells must map exactly to GSE81475 GEO records")
}
fetal_geo <- geo[geo_index, , drop = FALSE]
fetal_capture_device <- sub(
  "^((?:C1-)?[0-9]+)-.*$",
  "\\1",
  as.character(fetal_meta$Chip),
  perl = TRUE
)
fetal_capture_device <- ifelse(
  grepl("^C1-", fetal_capture_device),
  fetal_capture_device,
  paste0("C1-", fetal_capture_device)
)
if (
  length(unique(fetal_meta$Donor)) != 8L ||
    length(unique(fetal_capture_device)) != 12L ||
    any(vapply(
      split(as.character(fetal_meta$Donor), fetal_capture_device),
      function(value) length(unique(value)),
      integer(1)
    ) != 1L)
) {
  stop("Li fetal source must retain 8 donors across 12 donor-nested C1 chips")
}

adult_donor <- sub("DFC$", "", as.character(adult_meta$orig.ident))
adult_library_map <- c(
  HSB340 = "MS0314QR",
  HSB189 = "MS0314OP",
  HSB106 = "RQ7265KL"
)
adult_age_map <- c(HSB340 = "18.94 years", HSB189 = "36 years", HSB106 = "64 years")
adult_sex_map <- c(HSB340 = "M", HSB189 = "M", HSB106 = "M")
if (length(unique(adult_donor)) != 3L || any(!adult_donor %in% names(adult_library_map))) {
  stop("Li adult source must retain the three published DFC donors")
}

li_retained_scope <- paste(
  "two-arm retained source: adult raw UMI counts for the exact 17,093",
  "umi2/meta2 identities from 17,335 nuclei, plus 762 prenatal cells",
  "in the distributed PsychENCODE RData"
)
li_attrition_status <- paste(
  "242 adult nuclei are absent from the retained umi2/meta2 identities;",
  "the paper reports 1,195 prenatal cells but cell-level attrition to the",
  "762-cell distributed prenatal RData is not reconstructible"
)
li_diagnosis_status <- paste(
  "prenatal specimens and three adult DFC donors are recorded as",
  "unaffected/control in the paper and retained metadata"
)

make_fetal_metadata <- function() {
  cells <- colnames(fetal_counts)
  meta <- fetal_meta
  meta$Cells <- cells
  meta$Dataset <- "Li_et_al_2018"
  meta$Source_Li_Arm <- "fetal_scRNA"
  meta$Technology <- "Fluidigm C1"
  meta$Sequence <- "scRNA-seq"
  meta$Assay_Type <- "single-cell RNA sequencing"
  meta$Sequencing_Platform <- fetal_geo$Sample_instrument_model
  meta$Library_Chemistry <- paste(
    "Fluidigm C1 capture; SMARTer Ultra Low RNA cDNA;",
    "Illumina Nextera XT library"
  )
  meta$Original_Donor_ID <- as.character(meta$Donor)
  meta$Original_Specimen_ID <- paste(meta$Donor, meta$Region, sep = "::")
  meta$Original_Library_ID <- cells
  meta$Original_Source_Sample_ID <- meta$Original_Specimen_ID
  meta$Original_Source_Record_ID <- fetal_geo$Sample_geo_accession
  meta$Original_Technical_Batch_ID <- fetal_capture_device
  meta$Source_GEO_Record_ID <- fetal_geo$Sample_geo_accession
  meta$Source_Single_Cell_Library_ID <- cells
  meta$Source_Capture_Device_ID <- fetal_capture_device
  meta$Source_Cell_Type_Broad <- NA_character_
  meta$Source_Cell_Type_Fine <- as.character(meta$ctype)
  meta$CellType_raw <- meta$Source_Cell_Type_Fine
  meta$Source_Raw_Cell_Count <- NA_integer_
  meta$Source_Retained_Cell_Count <- 762L
  meta$Source_Excluded_Cell_Count <- NA_integer_
  meta$Source_Preprocessing_Retained_Object_Scope <- li_retained_scope
  meta$Source_Preprocessing_Attrition_Status <- li_attrition_status
  meta$Source_Preprocessing_Diagnosis_Audit_Status <- li_diagnosis_status
  meta$Sample <- meta$Original_Specimen_ID
  meta$Sample_ID <- meta$Original_Source_Record_ID
  meta$Brain_Region <- as.character(meta$Region)
  meta <- standardize_source_age_metadata(
    meta,
    source_age = as.character(meta$Age),
    source_basis = "postconceptional",
    source_unit = "post-conception weeks",
    source_reference = paste(
      "Li et al. 2018 original prenatal RData and Table S3,",
      "DOI 10.1126/science.aat7615"
    )
  )
  meta$Sex_Source_Raw <- as.character(meta$Sex)
  meta$Sex_Assignment_Method <- "source reported"
  meta$Sex <- standardize_source_sex_value(meta$Sex_Source_Raw)
  rownames(meta) <- cells
  meta
}

make_adult_metadata <- function() {
  cells <- colnames(adult_counts)
  meta <- adult_meta
  meta$Cells <- cells
  meta$Dataset <- "Li_et_al_2018"
  meta$Source_Li_Arm <- "adult_snRNA"
  meta$Technology <- "10x Genomics"
  meta$Sequence <- "snRNA-seq"
  meta$Assay_Type <- "single-nucleus RNA sequencing"
  meta$Sequencing_Platform <- NA_character_
  meta$Library_Chemistry <- "10x Genomics 3-prime gene expression; version not reported in retained source files"
  meta$Original_Donor_ID <- adult_donor
  meta$Original_Specimen_ID <- paste(adult_donor, "DFC", sep = "::")
  meta$Original_Library_ID <- unname(adult_library_map[adult_donor])
  meta$Original_Source_Sample_ID <- meta$Original_Specimen_ID
  meta$Original_Source_Record_ID <- meta$Original_Library_ID
  meta$Original_Technical_Batch_ID <- meta$Original_Library_ID
  meta$Source_10x_Library_ID <- meta$Original_Library_ID
  meta$Source_Cell_Type_Broad <- as.character(meta$ctype)
  meta$Source_Cell_Type_Fine <- as.character(meta$subtype)
  meta$CellType_raw <- meta$Source_Cell_Type_Fine
  meta$Source_Raw_Cell_Count <- 17335L
  meta$Source_Retained_Cell_Count <- 17093L
  meta$Source_Excluded_Cell_Count <- 242L
  meta$Source_Preprocessing_Retained_Object_Scope <- li_retained_scope
  meta$Source_Preprocessing_Attrition_Status <- li_attrition_status
  meta$Source_Preprocessing_Diagnosis_Audit_Status <- li_diagnosis_status
  meta$Sample <- meta$Original_Specimen_ID
  meta$Sample_ID <- meta$Original_Source_Record_ID
  meta$Brain_Region <- "Dorsolateral prefrontal cortex"
  meta$Region <- "DFC"
  meta <- standardize_source_age_metadata(
    meta,
    source_age = unname(adult_age_map[adult_donor]),
    source_unit = "years",
    source_reference = paste(
      "Li et al. 2018 original adult RData and Table S4,",
      "DOI 10.1126/science.aat7615"
    )
  )
  meta$Sex_Source_Raw <- unname(adult_sex_map[adult_donor])
  meta$Sex_Assignment_Method <- "source reported"
  meta$Sex <- standardize_source_sex_value(meta$Sex_Source_Raw)
  rownames(meta) <- cells
  meta
}

fetal_metadata <- make_fetal_metadata()
adult_metadata <- make_adult_metadata()
metadata_columns <- union(names(adult_metadata), names(fetal_metadata))
align_metadata <- function(meta) {
  for (column in setdiff(metadata_columns, names(meta))) {
    meta[[column]] <- NA
  }
  meta[, metadata_columns, drop = FALSE]
}
adult_metadata <- align_metadata(adult_metadata)
fetal_metadata <- align_metadata(fetal_metadata)
metadata <- rbind(adult_metadata, fetal_metadata)
metadata <- retain_source_metadata(
  metadata,
  c(
    "Cells", "Dataset", "Source_Li_Arm", "Technology", "Sequence",
    "Assay_Type", "Sequencing_Platform", "Library_Chemistry",
    "Original_Donor_ID", "Original_Specimen_ID", "Original_Library_ID",
    "Original_Source_Sample_ID", "Original_Source_Record_ID",
    "Original_Technical_Batch_ID", "Source_Cell_Type_Broad",
    "Source_Cell_Type_Fine", "CellType_raw", "Sample", "Sample_ID",
    "Source_Raw_Cell_Count", "Source_Retained_Cell_Count",
    "Source_Excluded_Cell_Count",
    "Source_Preprocessing_Retained_Object_Scope",
    "Source_Preprocessing_Attrition_Status",
    "Source_Preprocessing_Diagnosis_Audit_Status",
    "Brain_Region", "Region", "Age_Source_Raw", "Age_Source_Unit",
    "Age_Source_Basis", "Age_Source_Reference",
    "Age_Harmonization_Input", "Age_Conversion_Formula",
    "Age_Conversion_Confidence", "Age_Conversion_Applied", "Age",
    "Sex_Source_Raw", "Sex_Assignment_Method", "Sex"
  )
)

features <- union(rownames(adult_counts), rownames(fetal_counts))
adult_counts <- align_sparse_features(adult_counts, features)
fetal_counts <- align_sparse_features(fetal_counts, features)
input_cells <- c(colnames(adult_counts), colnames(fetal_counts))
metadata <- metadata[input_cells, , drop = FALSE]
if (!identical(metadata$Cells, input_cells)) {
  stop("Li et al. full metadata and matrix cell order differ")
}

thisutils::log_message("Creating the complete two-arm Li et al. 2018 object...")
object <- create_layered_processed_object(
  counts_layers = list(
    adult_snRNA = adult_counts,
    fetal_scRNA = fetal_counts
  ),
  metadata = metadata,
  dataset = "Li_et_al_2018",
  input_cells = input_cells,
  input_features = features
)
validate_processed_object(object, expected_cells = 17855L)
saveRDS(object, file.path(res_dir, "Li_et_al_2018_processed.rds"))
