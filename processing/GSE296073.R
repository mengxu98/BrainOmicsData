source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE296073/"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE296073/")

log_message("Start loading data...")

object <- readRDS(
  file.path(data_dir, "h_pre_peri_DY.rds")
)

object$Brain_Region <- recode(
  object$regionID,
  "cx" = "Cortex",
  "lv" = "Periventricular",
  "wh" = "Whole"
)
object$Age <- dplyr::case_when(
  object$age == "gw22" ~ "22 PCW",
  object$age == "gw23" ~ "23 PCW",
  object$age == "gw30" ~ "30 PCW",
  object$age == "pw02" ~ as.character(round(2 / 52, 2)),
  object$age == "pw03" ~ as.character(round(3 / 52, 2)),
  TRUE ~ NA_character_
)
object$Sex <- dplyr::recode(
  object$subjectID,
  "2023001" = "Male",
  "2023002" = "Female",
  "2024002" = "Male",
  "2014028" = "Female",
  "2018003" = "Male",
  "2018023" = "Female"
)

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
  file.path(res_dir, "GSE296073_processed.rds")
)
