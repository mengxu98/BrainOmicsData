source("functions/prepare_env.R")
source("functions/dataset_metadata.R")

res_dir <- check_dir(brainomics_data_path("processed/GSE186538/"))

thisutils::log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "GSE186538.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE186538"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$donor_ID
metadata$Sample_ID <- metadata$project_code
metadata$CellType_raw <- metadata$original_name
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$subregion
metadata <- standardize_source_age_metadata(
  metadata,
  source_age = as.character(metadata$donor_age),
  source_unit = "years",
  source_reference = "Franjic et al. 2022, DOI 10.1016/j.neuron.2021.10.036"
)
metadata$Sex_Source_Raw <- as.character(metadata$donor_gender)
metadata$Sex_Assignment_Method <- "source reported"
metadata$Sex <- standardize_source_sex_value(metadata$Sex_Source_Raw)

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
  file.path(res_dir, "GSE186538_processed.rds")
)
