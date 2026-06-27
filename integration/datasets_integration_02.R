library(Seurat)
library(thisutils)

data_dir <- "../../data/BrainOmicsData/integration/"

objects_list_file <- file.path(data_dir, "objects_list_processed.rds")
if (!file.exists(objects_list_file)) {
  log_message("Loading data...")
  objects_list <- readRDS(file.path(data_dir, "objects_list_raw.rds"))
  metadata <- readRDS(file.path(data_dir, "metadata_filtered.rds"))
  metadata_list <- split(metadata, metadata$Dataset)

  for (dataset in names(objects_list)) {
    log_message("Restricting dataset: {.val {dataset}}...")
    meta <- metadata_list[[dataset]]
    obj <- objects_list[[dataset]]
    obj <- obj[, meta$Cells]
    obj <- CreateSeuratObject(
      counts = GetAssayData(obj, layer = "counts"),
      meta.data = meta
    )
    obj <- JoinLayers(obj)
    objects_list[[dataset]] <- obj
    rm(obj, meta)
    gc()
  }
  log_message("Saving objects list...")
  saveRDS(objects_list, objects_list_file)
  log_message("Objects list saved to {.file {objects_list_file}}")
  rm(objects_list, metadata, metadata_list)
  gc()
}

objects_file <- file.path(data_dir, "objects_raw.rds")
if (!file.exists(objects_file)) {
  log_message("Loading objects list...")
  objects_list <- readRDS(file.path(data_dir, "objects_list_processed.rds"))

  log_message("Merging objects...")
  objects <- merge(
    objects_list[[1]],
    y = objects_list[-1]
  )
  rm(objects_list)
  gc()

  objects[["percent.mt"]] <- PercentageFeatureSet(
    objects,
    pattern = "^MT-"
  )
  log_message("Saving objects...")
  saveRDS(objects, objects_file)
  log_message("Objects saved to {.file {objects_file}}")
}

objects_file <- file.path(data_dir, "objects_filtered.rds")
if (!file.exists(objects_file)) {
  objects <- readRDS(file.path(data_dir, "objects_raw.rds"))
  remove_patterns <- c(
    "^ERCC", "^RPLP", "^RPSL", "^MT-", "^mt-",
    "^LOC", "^LINC", "^RP[0-9]", "^AC[0-9]", "MALAT1", "^HB[^(P)]",
    "^AL[0-9]", "^AP[0-9]", "^CT[0-9]", "^CH[0-9]", "^FAM[0-9]", "orf"
  )
  keep <- !grepl(
    paste(
      remove_patterns,
      collapse = "|"
    ), rownames(objects),
    ignore.case = TRUE
  )
  objects <- objects[keep, ]
  saveRDS(objects, objects_file)
  log_message("Objects saved to {.file {objects_file}}")
}
