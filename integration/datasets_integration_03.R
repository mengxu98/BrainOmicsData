library(Seurat)
library(thisutils)
library(future)
library(harmony)

options(future.globals.maxSize = 8 * 1024^3)
plan(sequential)

dims <- 1:50
n_features <- 3000

res_dir <- "../../data/BrainOmicsData/integration/"

objects_file <- file.path(
  res_dir, "objects_integrated.rds"
)
if (!file.exists(objects_file)) {
  log_message("Number of parallel workers: {.val {nbrOfWorkers()}}")

  log_message("Loading objects")
  objects <- readRDS(file.path(res_dir, "objects_filtered.rds"))

  log_message("Running preprocessing")

  objects <- NormalizeData(objects)
  objects <- FindVariableFeatures(
    objects,
    nfeatures = n_features
  )
  objects <- ScaleData(objects)
  objects <- RunPCA(
    objects,
    features = VariableFeatures(objects)
  )

  objects <- FindNeighbors(
    objects,
    dims = dims,
    reduction = "pca"
  )
  objects <- FindClusters(
    objects,
    resolution = 1
  )
  objects <- RunUMAP(
    objects,
    dims = dims,
    reduction = "pca",
    reduction.name = "umap.unintegrated"
  )

  log_message("Integrating layers using RPCA")
  objects <- IntegrateLayers(
    object = objects,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca"
  )
  log_message("Running UMAP using RPCA")
  objects <- RunUMAP(
    objects,
    reduction = "integrated.rpca",
    dims = dims,
    reduction.name = "umap.rpca"
  )

  log_message("Integrating layers using Harmony")
  objects <- IntegrateLayers(
    object = objects,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.harmony"
  )
  log_message("Running UMAP using Harmony")
  objects <- RunUMAP(
    objects,
    reduction = "integrated.harmony",
    dims = dims,
    reduction.name = "umap.harmony"
  )

  pca_raw <- objects@reductions$pca@cell.embeddings
  pca_rpca <- objects@reductions$integrated.rpca@cell.embeddings
  pca_harmony <- objects@reductions$integrated.harmony@cell.embeddings
  umap_raw <- objects@reductions$umap.unintegrated@cell.embeddings
  umap_rpca <- objects@reductions$umap.rpca@cell.embeddings
  umap_harmony <- objects@reductions$umap.harmony@cell.embeddings
  meta_data <- objects@meta.data
  lisi_data <- list(
    pca_raw,
    pca_rpca,
    pca_harmony,
    umap_raw,
    umap_rpca,
    umap_harmony,
    meta_data
  )
  names(lisi_data) <- c(
    "pca_raw", "pca_rpca", "pca_harmony",
    "umap_raw", "umap_rpca", "umap_harmony",
    "meta_data"
  )
  saveRDS(
    lisi_data,
    file.path(res_dir, "lisi_data.rds")
  )

  log_message("Saving objects")
  saveRDS(
    objects,
    objects_file
  )
}
