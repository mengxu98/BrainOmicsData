source("functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE186538/")

log_message("Start loading data...")
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
metadata$Age <- metadata$donor_age
metadata$Age <- gsub("yr", "", metadata$Age)
metadata$Sex <- metadata$donor_gender
metadata$Sex <- ifelse(metadata$Sex == "F", "Female", "Male")

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- GetAssayData(object, layer = "counts")
counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE186538_processed.rds")
)
