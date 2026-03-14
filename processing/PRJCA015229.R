source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/PRJCA015229"
res_dir <- check_dir("../../data/BrainOmicsData/processed/PRJCA015229/")

log_message("Start loading data...")

data_path <- read.csv(
  file.path(data_dir, "PRJCA015229/data.detail.csv")
)
file_object_list <- file.path(res_dir, "object_list(6samples)_snRNA-seq.rds")
if (!file.exists(file_object_list)) {
  data_list <- list()
  data_path_multiomes <- subset(data_path, data.type == "multiomes")
  data_path_multiomes <- subset(data_path_multiomes, species == "Human")
  data_path_rna <- subset(data_path, data.type == "snRNA")
  data_path_rna <- subset(data_path_rna, species == "Human")

  for (i in seq_along(data_path_multiomes$objectID)) {
    data_list[[i]] <- Read10X(
      file.path(
        data_dir,
        gsub("data/", "", data_path_multiomes$pathway[i]),
        "filtered_feature_bc_matrix"
      ),
      gene.column = 2
    )
    data_list[[i]] <- data_list[[i]][[1]]
  }
  for (i in seq_along(data_path_rna$objectID)) {
    data_list[[i + length(data_path_multiomes$objectID)]] <- Read10X(
      file.path(
        data_dir,
        gsub("data/", "", data_path_rna$pathway[i]),
        "filtered_feature_bc_matrix"
      ),
      gene.column = 2
    )
  }

  human_id <- purrr::map(
    data_list, function(x) unique(rownames(x))
  ) |>
    purrr::list_c() |>
    unique()

  for (i in seq_along(data_list)) {
    data_list[[i]] <- data_list[[i]][human_id, ]
  }
  names(data_list) <- c(data_path_multiomes$objectID, data_path_rna$objectID)

  object_list <- list()
  for (i in seq_along(data_list)) {
    object_list[[i]] <- CreateSeuratObject(
      counts = data_list[[i]],
      project = names(data_list)[i]
    )
    object_list[[i]][["homo"]] <- CreateAssayObject(
      counts = data_list[[i]]
    )
  }
  names(object_list) <- names(data_list)
  rm(data_list)
  gc()

  saveRDS(
    object_list,
    file.path(
      file_object_list
    )
  )

  for (i in seq_along(object_list)) {
    object_list[[i]] <- scop::RunDoubletCalling(
      object_list[[i]],
      db_method = "Scrublet"
    )
    object_list[[i]]@meta.data[["DoubletScores"]] <- object_list[[i]]@meta.data[["db.Scrublet_score"]]
    object_list[[i]]@meta.data[["PredictedDoublets"]] <- object_list[[i]]@meta.data[["db.Scrublet_class"]]
    object_list[[i]]@meta.data[["DoubletScores"]] <- unlist(object_list[[i]]@meta.data[["DoubletScores"]])
    object_list[[i]]@meta.data[["PredictedDoublets"]] <- unlist(object_list[[i]]@meta.data[["PredictedDoublets"]])
  }

  for (i in seq_along(object_list)) {
    object_list[[i]][["percent_mito"]] <- PercentageFeatureSet(
      object_list[[i]],
      pattern = "^MT-"
    )
  }
  object_list <- lapply(
    X = object_list, FUN = function(x) {
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(
        x,
        selection.method = "vst", nfeatures = 2000, verbose = FALSE
      )
      x <- ScaleData(x, verbose = FALSE)
      x <- RunPCA(x, npcs = 30, verbose = FALSE)
      x <- RunUMAP(x, reduction = "pca", dims = 1:30, verbose = FALSE)
      x <- FindNeighbors(x, reduction = "pca", dims = 1:30, verbose = FALSE)
      x <- FindClusters(x, resolution = 0.5, verbose = FALSE)
    }
  )

  saveRDS(
    object_list,
    file.path(
      file_object_list
    )
  )
} else {
  object_list <- readRDS(
    file.path(
      file_object_list
    )
  )
}

object <- merge(object_list[[1]], object_list[2:6])
object <- JoinLayers(object)

object <- subset(
  object,
  PredictedDoublets == "singlet"
)

cell_ids <- gsub("_[0-9]+$", "", colnames(object))
orig_ident <- object$orig.ident
new_colnames <- paste0(cell_ids, "_", orig_ident)
colnames(object) <- new_colnames

# https://www.cell.com/cms/10.1016/j.xgen.2024.100703/attachment/29b2bce9-946b-43b2-8f12-346dbaea1b6f/mmc6.xlsx
metadata_snmultiome <- read.csv(
  file.path(data_dir, "metadata_snMultiome.csv")
)
metadata_snrnaseq <- read.csv(
  file.path(data_dir, "metadata_snRNA-seq.csv")
)
common_colnames <- intersect(
  colnames(metadata_snmultiome), colnames(metadata_snrnaseq)
)
metadata_snrnaseq$Barcode <- paste0(
  metadata_snrnaseq$Barcode, "_", metadata_snrnaseq$Sample
)
metadata <- rbind(
  metadata_snmultiome[, common_colnames],
  metadata_snrnaseq[, common_colnames]
)

common_cells <- intersect(new_colnames, metadata$Barcode)
log_message(
  "Common cells between object and metadata: {.val {length(common_cells)}}"
)

object <- object[, common_cells]

metadata_filtered <- metadata[metadata$Barcode %in% common_cells, ]
index <- match(colnames(object), metadata_filtered$Barcode)
metadata_filtered <- metadata_filtered[index, ]
rownames(metadata_filtered) <- metadata_filtered$Barcode

metadata <- metadata_filtered
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "PRJCA015229"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample_ID <- metadata$Sample
metadata$CellType_raw <- metadata$BigCellType
metadata$Brain_Region <- "Anterior cingulate cortex "
metadata$Region <- "Anterior cingulate cortex "

sample_info <- data.frame(
  Sample = c(
    "HM2013017", "HM20200905", "HM20200927", "HM20201129", "HM20201213", "HM20201222"
  ),
  Age = c(
    "58", "44", "52", "69", "47", "66"
  ),
  Sex = c(
    "Male", "Male", "Male", "Male", "Male", "Male"
  ),
  stringsAsFactors = FALSE
)

metadata <- merge(
  metadata, sample_info,
  by.x = "Sample", by.y = "Sample", all.x = TRUE
)
rownames(metadata) <- metadata$Cells

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

saveRDS(
  object,
  file.path(
    res_dir, "PRJCA015229_processed.rds"
  )
)
