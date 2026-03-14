source("functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainOmicsData/processed/Nowakowski_et_al_2017/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "Nowakowski_et_al_2017.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "Nowakowski_et_al_2017"
metadata$Technology <- metadata$seq_tech
metadata$Sequence <- metadata$seq_method
metadata$Sample <- metadata$donor_ID
metadata$Sample_ID <- metadata$orig.ident
metadata$CellType_raw <- metadata$original_name
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$region
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
  file.path(res_dir, "Nowakowski_et_al_2017_processed.rds")
)
