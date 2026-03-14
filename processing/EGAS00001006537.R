source("functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainOmicsData/processed/EGAS00001006537/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "EGAS00001006537_raw.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "EGAS00001006537"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$sampleID
metadata$Sample_ID <- metadata$sampleID
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
  file.path(res_dir, "EGAS00001006537_processed.rds")
)
