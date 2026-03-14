source("functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE104276/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "GSE104276.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE104276"
metadata$Technology <- "STRT-seq"
metadata$Sequence <- "scRNA-seq"
metadata$Sample <- metadata$sample
metadata$Sample_ID <- metadata$cell_ID
metadata$CellType_raw <- metadata$cell_types
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$subregion
metadata$Age <- metadata$donor_age
metadata$Age <- gsub("w", " PCW", metadata$Age)
metadata <- metadata[metadata$Age != "Unclassified", ]
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
  file.path(res_dir, "GSE104276_processed.rds")
)
