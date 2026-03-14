source("functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainOmicsData/processed/EGAD00001006049/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "EGAD00001006049_raw.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "EGAD00001006049"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- metadata$Library
metadata$Sample <- metadata$donorID
metadata$Sample_ID <- metadata$donorID
metadata$CellType_raw <- metadata$Cell.Type
metadata$Brain_Region <- metadata$Regions
metadata$Region <- metadata$Regions
metadata$Age_new <- sapply(
  seq_len(nrow(metadata)), function(i) {
    stage <- metadata$Stage[i]
    age <- metadata$Age[i]

    if (grepl("Fetal", stage)) {
      paste0(round(age / 7, 1), " PCW")
    } else {
      round(age / 365, 0)
    }
  }
)
metadata$Age <- metadata$Age_new
metadata$Sex <- ifelse(metadata$Sex_impute == "F", "Female", "Male")

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
  file.path(res_dir, "EGAD00001006049_processed.rds")
)
