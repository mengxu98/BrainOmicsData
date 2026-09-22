source("functions/prepare_env.R")
source("functions/dataset_metadata.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE104276/rawData"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE104276/")

thisutils::log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "GSE104276.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE104276"
metadata$Technology <- "modified Smart-seq2 with UMI barcodes"
metadata$Sequence <- "scRNA-seq"
metadata$Assay_Type <- "single-cell RNA-seq"
metadata$Sequencing_Platform <- "Illumina HiSeq 4000"
metadata$Library_Chemistry <- "modified Smart-seq2 with UMI barcodes"

source_cell_id <- as.character(metadata$cell_ID)
source_donor_id <- sub(
  "^(GW[0-9]+_PFC[0-9]+).*$",
  "\\1",
  source_cell_id
)
if (
  any(!grepl("^GW[0-9]+_PFC[0-9]+$", source_donor_id)) ||
    length(unique(source_donor_id)) != 12L
) {
  stop("GSE104276 must retain the 12 publication-defined PFC donors")
}

source_tar <- file.path(data_dir, "GSE104276_RAW.tar")
if (!file.exists(source_tar)) {
  stop("GSE104276 GEO raw archive is missing: ", source_tar)
}
source_members <- utils::untar(source_tar, list = TRUE)
if (length(source_members) != 39L) {
  stop("GSE104276 GEO raw archive must contain 39 GSM records")
}
extract_dir <- tempfile("gse104276-geo-")
dir.create(extract_dir)
on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
utils::untar(source_tar, exdir = extract_dir)
source_record_map <- do.call(
  rbind,
  lapply(source_members, function(member) {
    connection <- gzfile(file.path(extract_dir, member), "rt")
    header <- readLines(connection, n = 1L)
    close(connection)
    cells <- strsplit(header, "\t", fixed = TRUE)[[1L]][-1L]
    data.frame(
      cell_ID = cells,
      Source_GEO_Record_ID = sub("_.*", "", basename(member)),
      Source_GEO_Record_File = basename(member),
      stringsAsFactors = FALSE
    )
  })
)
if (
  anyDuplicated(source_record_map$cell_ID) ||
    length(unique(source_record_map$Source_GEO_Record_ID)) != 39L
) {
  stop("GSE104276 GEO cell-to-record mapping is malformed")
}

# Ten retained combined-matrix labels occur exactly at the boundary between
# consecutive GEO records and are omitted from the individual-file headers.
# Their record identity follows the published sample/barcode workbook and the
# immediately following record in each boundary pair.
boundary_record_map <- c(
  GW10_PFC2_sc49 = "GSM2884064",
  GW19_PFC1_sc41 = "GSM2884078",
  GW19_PFC1_sc81 = "GSM2884079",
  GW23_PFC1_sc49 = "GSM2884081",
  GW23_PFC1_sc97 = "GSM2884082",
  GW23_PFC2_SF1_F23_sc26 = "GSM2970391",
  GW23_PFC2_SF2_F25_sc26 = "GSM2970392",
  GW23_PFC2_PAO1_F24_sc1 = "GSM2970393",
  GW23_PFC2_RMF2_F16_sc26 = "GSM2970394",
  GW23_PFC2_PAO3_F15_sc27 = "GSM2970395"
)
source_record <- source_record_map$Source_GEO_Record_ID[
  match(source_cell_id, source_record_map$cell_ID)
]
missing_record <- is.na(source_record)
source_record[missing_record] <- unname(
  boundary_record_map[source_cell_id[missing_record]]
)
if (
  anyNA(source_record) ||
    length(unique(source_record)) != 39L ||
    any(vapply(
      split(source_donor_id, source_record),
      function(value) length(unique(value)) != 1L,
      logical(1)
    ))
) {
  stop("GSE104276 retained cells do not map to 39 donor-specific GSM records")
}

metadata$Sample <- source_donor_id
metadata$Sample_ID <- source_record
metadata$Original_Source_Sample_ID <- source_donor_id
metadata$Original_Source_Record_ID <- source_record
metadata$Source_GEO_Record_ID <- source_record
metadata$Source_Single_Cell_Library_ID <- source_cell_id
metadata$Source_Published_Cell_Count <- 2394L
metadata$Source_Retained_Cell_Count <- nrow(metadata)
metadata$Source_Cell_Type_Broad_Original <- as.character(
  metadata$cell_types
)
metadata$Source_Secondary_Cell_Type_Label <- as.character(
  metadata$cell_type
)
metadata$Source_Secondary_Cell_Type_Semantics <- paste(
  "secondary reference label retained from the input object;",
  "not substituted for a missing original-study cell_types label"
)
metadata$CellType_raw <- metadata$cell_types
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$subregion
source_age_original <- as.character(metadata$donor_age)
source_age_from_cell_id <- paste0(
  as.integer(sub("^GW([0-9]+)_.*$", "\\1", source_cell_id)),
  "w"
)
classified_age <- !is.na(source_age_original) &
  source_age_original != "" &
  source_age_original != "Unclassified"
if (any(source_age_original[classified_age] !=
  source_age_from_cell_id[classified_age])) {
  stop("GSE104276 donor-age labels disagree with the source cell IDs")
}
source_age <- source_age_original
age_reconstructed <- is.na(source_age) |
  source_age == "" |
  source_age == "Unclassified"
source_age[age_reconstructed] <- source_age_from_cell_id[age_reconstructed]
metadata$Source_Donor_Age_Original <- source_age_original
metadata$Source_Age_From_Cell_ID <- source_age_from_cell_id
metadata$Source_Age_Derivation_Rule <- ifelse(
  age_reconstructed,
  paste(
    "replace the non-age label Unclassified with the gestational week",
    "encoded by the publication cell-ID GW## prefix"
  ),
  "retain donor_age after exact agreement with the cell-ID GW## prefix"
)
gestational <- !is.na(source_age) & grepl("w", source_age, ignore.case = TRUE)
metadata <- standardize_source_age_metadata(
  metadata,
  source_age = source_age,
  source_basis = ifelse(gestational, "gestational", NA_character_),
  source_unit = ifelse(
    gestational,
    "gestational weeks",
    "source developmental-age label"
  ),
  source_reference = paste(
    "Zhong et al. 2018, DOI 10.1038/nature25980 and GEO GSE104276;",
    "the source design explicitly reports gestational weeks 8-26"
  )
)
metadata$Sex_Source_Raw <- as.character(metadata$donor_gender)
metadata$Sex_Assignment_Method <- "source reported"
metadata$Sex <- standardize_source_sex_value(metadata$Sex_Source_Raw)

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "Original_Source_Sample_ID",
  "Original_Source_Record_ID", "Source_GEO_Record_ID",
  "Source_Single_Cell_Library_ID", "Source_Published_Cell_Count",
  "Source_Retained_Cell_Count", "Assay_Type", "Sequencing_Platform",
  "Library_Chemistry", "Source_Cell_Type_Broad_Original",
  "Source_Secondary_Cell_Type_Label",
  "Source_Secondary_Cell_Type_Semantics", "CellType_raw",
  "Brain_Region", "Region", "Source_Donor_Age_Original",
  "Source_Age_From_Cell_ID", "Source_Age_Derivation_Rule",
  "Age_Source_Raw", "Age_Source_Unit", "Age_Source_Basis",
  "Age_Source_Reference", "Age_Harmonization_Input",
  "Age_Conversion_Formula", "Age_Conversion_Confidence",
  "Age_Conversion_Applied", "Age", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)

counts <- GetAssayData(object, layer = "counts")
counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

thisutils::log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE104276_processed.rds")
)
