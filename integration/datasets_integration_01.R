source("functions/prepare_env.R")
source("functions/metadata_schema.R")

data_dir <- "../../data/BrainOmicsData/processed/"
res_dir <- check_dir(
  "../../data/BrainOmicsData/integration/"
)

log_message("Start loading data...")

datasets <- c(
  "AllenM1", "EGAD00001006049", "EGAS00001006537",
  "GSE103723", "GSE104276", "GSE144136", "GSE168408",
  "GSE178175", "GSE186538", "GSE199762", "GSE202210",
  "GSE204683", "GSE207334", "GSE212606", "GSE217511",
  "GSE261983", "GSE296073", "GSE67835", "GSE81475", "GSE97942",
  "HYPOMAP", "Li_et_al_2018",
  "Ma_et_al_2022", "Nowakowski_et_al_2017",
  "PRJCA015229", "ROSMAP", "SomaMut"
)

objects_list_file <- file.path(res_dir, "objects_list_raw.rds")
if (!file.exists(objects_list_file)) {
  objects_list <- list()
  for (dataset in datasets) {
    log_message("Loading dataset: {.val {dataset}}...")
    file_path <- file.path(
      data_dir, dataset, paste0(dataset, "_processed.rds")
    )
    object <- readRDS(file_path)
    object <- JoinLayers(object)
    objects_list[[dataset]] <- object
    rm(object)
    gc()
  }
  names(objects_list) <- datasets

  metadata <- purrr::map_dfr(
    objects_list,
    function(x) {
      meta <- x@meta.data
      for (col in colnames(meta)) {
        meta[[col]] <- as.character(meta[[col]])
      }
      meta
    }
  )

  metadata$Unit <- ifelse(grepl("PCW", metadata$Age), "PCW", "Years")
  get_num <- function(x) {
    as.numeric(sub("([0-9.]+).*", "\\1", x))
  }
  metadata$Age_num <- NA_real_
  pcw_idx <- metadata$Unit == "PCW"
  metadata$Age_num[pcw_idx] <- get_num(metadata$Age[pcw_idx])
  yr_idx <- metadata$Unit == "Years"
  metadata$Age_num[yr_idx] <- as.numeric(metadata$Age[yr_idx])
  metadata$DevelopmentStage <- dplyr::case_when(
    metadata$Age_num >= 4 & metadata$Age_num < 8 & metadata$Unit == "PCW" ~ "Embryonic (4-8 PCW)",
    metadata$Age_num >= 8 & metadata$Age_num < 10 & metadata$Unit == "PCW" ~ "Early fetal (8-10 PCW)",
    metadata$Age_num >= 10 & metadata$Age_num < 13 & metadata$Unit == "PCW" ~ "Early fetal (10-13 PCW)",
    metadata$Age_num >= 13 & metadata$Age_num < 16 & metadata$Unit == "PCW" ~ "Early mid-fetal (13-16 PCW)",
    metadata$Age_num >= 16 & metadata$Age_num < 19 & metadata$Unit == "PCW" ~ "Early mid-fetal (16-19 PCW)",
    metadata$Age_num >= 19 & metadata$Age_num < 24 & metadata$Unit == "PCW" ~ "Late mid-fetal (19-24 PCW)",
    metadata$Age_num >= 24 & metadata$Age_num < 39 & metadata$Unit == "PCW" ~ "Late fetal (24-38 PCW)",
    metadata$Age_num >= 0 & metadata$Age_num < 0.5 & metadata$Unit == "Years" ~ "Neonatal and early infancy (0-0.5Y)",
    metadata$Age_num >= 0.5 & metadata$Age_num < 1 & metadata$Unit == "Years" ~ "Late infancy (0.5-1Y)",
    metadata$Age_num >= 1 & metadata$Age_num < 6 & metadata$Unit == "Years" ~ "Early childhood (1-6Y)",
    metadata$Age_num >= 6 & metadata$Age_num < 12 & metadata$Unit == "Years" ~ "Middle and late childhood (6-12Y)",
    metadata$Age_num >= 12 & metadata$Age_num < 20 & metadata$Unit == "Years" ~ "Adolescence (12-20Y)",
    metadata$Age_num >= 20 & metadata$Age_num < 40 & metadata$Unit == "Years" ~ "Young adulthood (20-40Y)",
    metadata$Age_num >= 40 & metadata$Age_num < 60 & metadata$Unit == "Years" ~ "Middle adulthood (40-60Y)",
    metadata$Age_num >= 60 & metadata$Unit == "Years" ~ "Late adulthood (60+Y)",
    TRUE ~ NA_character_
  )

  metadata$Stage <- dplyr::case_when(
    metadata$Age_num >= 4 & metadata$Age_num < 8 & metadata$Unit == "PCW" ~ "S1",
    metadata$Age_num >= 8 & metadata$Age_num < 10 & metadata$Unit == "PCW" ~ "S2",
    metadata$Age_num >= 10 & metadata$Age_num < 13 & metadata$Unit == "PCW" ~ "S3",
    metadata$Age_num >= 13 & metadata$Age_num < 16 & metadata$Unit == "PCW" ~ "S4",
    metadata$Age_num >= 16 & metadata$Age_num < 19 & metadata$Unit == "PCW" ~ "S5",
    metadata$Age_num >= 19 & metadata$Age_num < 24 & metadata$Unit == "PCW" ~ "S6",
    metadata$Age_num >= 24 & metadata$Age_num < 39 & metadata$Unit == "PCW" ~ "S7",
    metadata$Age_num >= 0 & metadata$Age_num < 0.5 & metadata$Unit == "Years" ~ "S8",
    metadata$Age_num >= 0.5 & metadata$Age_num < 1 & metadata$Unit == "Years" ~ "S9",
    metadata$Age_num >= 1 & metadata$Age_num < 6 & metadata$Unit == "Years" ~ "S10",
    metadata$Age_num >= 6 & metadata$Age_num < 12 & metadata$Unit == "Years" ~ "S11",
    metadata$Age_num >= 12 & metadata$Age_num < 20 & metadata$Unit == "Years" ~ "S12",
    metadata$Age_num >= 20 & metadata$Age_num < 40 & metadata$Unit == "Years" ~ "S13",
    metadata$Age_num >= 40 & metadata$Age_num < 60 & metadata$Unit == "Years" ~ "S14",
    metadata$Age_num >= 60 & metadata$Unit == "Years" ~ "S15",
    TRUE ~ NA_character_
  )
  metadata <- add_metadata_schema(metadata)
  metadata <- metadata[metadata$BrainRegion != "Choroid", ]

  colnames_order <- c(
    "Cells", "Dataset", "Technology", "Sequence", "Sample", "Sample_ID",
    "CellType_raw", "BrainRegion", "Stage", "DevelopmentStage",
    "AgeIntervalID", "AgeInterval", "AgeRange", "Age", "Sex",
    "sequencing_modality_raw", "sequencing_modality_standardized",
    "sequencing_technology", "age_raw", "age_value", "age_unit", "age_sort",
    "sex_raw", "sex_standardized"
  )
  metadata$CellType_raw[is.na(metadata$CellType_raw)] <- "Unknown"

  metadata <- metadata[, colnames_order]
  metadata <- na.omit(metadata)
  saveRDS(metadata, file.path(res_dir, "metadata_filtered.rds"))
  saveRDS(objects_list, objects_list_file)
}
