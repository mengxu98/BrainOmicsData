# paper: https://doi.org/10.1038/s12276-024-01328-6
# data: https://zenodo.org/records/10939707


source("functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainOmicsData/processed/BTSatlas/")

log_message("Start loading data...")

if (!file.exists(file.path(res_dir, "BTS_atlas_raw.rds"))) {
  PrepareEnv()
  h5ad_path <- file.path(res_dir, "BTS_atlas_compatible.h5ad")

  object <- scop::h5ad_to_srt(h5ad_path, verbose = TRUE)
  gc()
  saveRDS(object, file.path(res_dir, "BTS_atlas_raw.rds"))
} else {
  object <- readRDS(file.path(res_dir, "BTS_atlas_raw.rds"))
}

dataset_info <- data.frame(
  Dataset = c(
    "AllenM1", "Braun", "Cameron", "Hardwick",
    "Herring", "Morabito", "Nagy", "Zhu"
  ),
  Dataset_ID = c(
    "AllenM1", "EGAD00001006049", "EGAS00001006537",
    "GSE178175", "GSE168408", "GSE174367",
    "GSE144136", "GSE202210"
  ),
  stringsAsFactors = FALSE
)

object$Dataset_ID <- dataset_info$Dataset_ID[match(
  object$Dataset, dataset_info$Dataset
)]

region_map <- c(
  # Prefrontal / Frontal Cortex
  "BA8" = "Prefrontal cortex",
  "BA9" = "Prefrontal cortex",
  "BA10" = "Prefrontal cortex",
  "BA46" = "Prefrontal cortex",
  "Frontal Cortex" = "Prefrontal cortex",
  # Motor cortex
  "M1" = "Primary motor cortex",
  # Temporal / Entorhinal
  "Cortex temporal" = "Temporal cortex",
  "Cortex entorhinal" = "Entorhinal cortex",
  # General cortex
  "Cortex" = "Cerebral cortex",
  "Brain" = "Cerebral cortex",
  "Forebrain" = "Cerebral cortex",
  "Subcortex" = "Subcortical region",
  # Basal ganglia
  "Caudate+Putamen" = "Basal ganglia",
  "Ganglionic Eminence" = "Ganglionic eminence",
  # Diencephalon / Thalamus / Hypothalamus
  "Diencephalon" = "Diencephalon",
  "Thalamus" = "Thalamus",
  "Hypothalamus" = "Hypothalamus",
  # Hindbrain & related
  "Hindbrain" = "Hindbrain",
  "Cerebellum" = "Cerebellum",
  "Pons" = "Pons",
  "Medulla" = "Medulla",
  # Midbrain
  "Midbrain" = "Midbrain",
  "Midbrain ventral" = "Midbrain ventral",
  # Limbic / Hippocampus
  "Hippocampus" = "Hippocampus"
)
brain_region_vec <- as.character(object$Brain.Region)
mapped_regions <- region_map[brain_region_vec]
names(mapped_regions) <- colnames(object)
object$Regions <- mapped_regions
na_indices <- is.na(object$Regions)
object$Regions[na_indices] <- brain_region_vec[na_indices]

cells_to_keep <- object$Dataset_ID != "GSE174367" # Alzheimer disease

object <- object[, cells_to_keep]


for (dataset in unique(object$Dataset_ID)) {
  log_message("Processing {.val {dataset}}...")
  dir_sub <- check_dir(
    file.path("../../data/BrainOmicsData/processed", dataset)
  )
  file_sub <- file.path(dir_sub, paste0(dataset, "_raw.rds"))
  if (file.exists(file_sub)) {
    log_message("File already exists: {.val {file_sub}}")
    next
  }
  cells_to_keep <- object$Dataset_ID == dataset
  object_dataset <- object[, cells_to_keep]
  saveRDS(object_dataset, file_sub)
  log_message("Saved {.val {dataset}}")
}
