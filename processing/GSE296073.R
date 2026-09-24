source("functions/prepare_env.R")
source("functions/dataset_metadata.R")

data_dir <- brainomics_data_path("raw/GSE296073/")
res_dir <- check_dir(brainomics_data_path("processed/GSE296073/"))

thisutils::log_message("Start loading data...")

object <- readRDS(
  file.path(data_dir, "h_pre_peri_DY.rds")
)

region_code <- as.character(object$regionID)
region_map <- c(
  cx = "Cortex",
  lv = "Periventricular",
  wh = "Whole"
)
if (any(!is.na(region_code) & !region_code %in% names(region_map))) {
  stop("GSE296073 contains an unsupported source region code")
}
object$Brain_Region <- unname(region_map[region_code])
age_code <- as.character(object$age)
gestational <- grepl("^gw[0-9]+$", age_code, ignore.case = TRUE)
postnatal_week <- grepl("^pw[0-9]+$", age_code, ignore.case = TRUE)
age_number <- suppressWarnings(as.numeric(sub("^[gp]w", "", age_code)))
if (any(!is.na(age_code) & !(gestational | postnatal_week))) {
  stop("GSE296073 contains an unsupported source age code")
}
canonical_age_value <- ifelse(
  gestational,
  age_number - 2,
  age_number / 52
)
canonical_age_unit <- ifelse(gestational, "PCW", "years")
object$Age_Source_Raw <- age_code
object$Age_Source_Unit <- ifelse(
  gestational,
  "gestational weeks",
  "postnatal weeks"
)
object$Age_Source_Basis <- ifelse(
  gestational,
  "gestational age",
  "postnatal age"
)
object$Age_Source_Reference <- paste(
  "Yang et al. 2025, DOI 10.1038/s41586-025-09362-8 and GEO GSE296073;",
  "source age codes use gw for gestational week and pw for postnatal week"
)
object$Age_Harmonization_Input <- paste(
  format_canonical_age_number(canonical_age_value),
  canonical_age_unit
)
object$Age_Conversion_Formula <- ifelse(
  gestational,
  "source gw - 2 weeks = PCW",
  "source postnatal weeks / 52 = years"
)
object$Age_Conversion_Confidence <-
  "high: explicit source age code and publication-linked metadata"
object$Age_Conversion_Applied <- TRUE
object$Age <- object$Age_Harmonization_Input
object$Sex_Derivation_Input <- as.character(object$subjectID)
object$Sex <- dplyr::recode(
  sub("^h", "", as.character(object$subjectID)),
  "2023001" = "Male",
  "2023002" = "Female",
  "2024002" = "Male",
  "2014028" = "Female",
  "2018003" = "Male",
  "2018023" = "Female",
  .default = NA_character_
)
if (anyNA(object$Sex)) {
  stop("GSE296073 subject-to-sex mapping is incomplete")
}
object$Sex_Source_Raw <- as.character(object$Sex)
object$Sex_Source_Standardized <- as.character(object$Sex)
object$Sex_Source_Column <-
  "source publication sample table keyed by subjectID"
object$Sex_Assignment_Method <-
  "curated from source publication sample table; not inferred"

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE296073"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$subjectID
metadata$Sample_ID <- metadata$libraryID
metadata$Region <- metadata$Brain_Region
metadata$CellType_raw <- metadata$cluster1

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region",
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
  file.path(res_dir, "GSE296073_processed.rds")
)
