source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/PRJCA015229"
res_dir <- check_dir("../../data/BrainOmicsData/processed/PRJCA015229/")

thisutils::log_message("Start loading data...")

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

cell_ids <- gsub("_[0-9]+$", "", colnames(object))
orig_ident <- object$orig.ident
new_colnames <- paste0(cell_ids, "_", orig_ident)
colnames(object) <- new_colnames

object_metadata <- object[[]]
if (!all(c(
  "db.Scrublet_score", "db.Scrublet_class", "PredictedDoublets"
) %in% names(object_metadata))) {
  stop("PRJCA015229 source Scrublet fields are missing")
}

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
thisutils::log_message(
  "Common cells between object and metadata: {.val {length(common_cells)}}"
)

object <- object[, common_cells]

metadata_filtered <- metadata[metadata$Barcode %in% common_cells, ]
index <- match(colnames(object), metadata_filtered$Barcode)
metadata_filtered <- metadata_filtered[index, ]
rownames(metadata_filtered) <- metadata_filtered$Barcode

metadata <- metadata_filtered
object_index <- match(metadata$Barcode, rownames(object_metadata))
if (anyNA(object_index)) {
  stop("PRJCA015229 Scrublet metadata does not cover source-annotated nuclei")
}
metadata$Source_Scrublet_Score <- as.numeric(
  object_metadata$db.Scrublet_score[object_index]
)
metadata$Source_Scrublet_Class <- as.character(
  object_metadata$db.Scrublet_class[object_index]
)
metadata$Source_Predicted_Doublet <- as.character(
  object_metadata$PredictedDoublets[object_index]
)
if (nrow(metadata) != 40748L ||
  sum(metadata$Source_Scrublet_Class == "doublet") != 576L ||
  sum(metadata$Source_Scrublet_Class == "singlet") != 40172L ||
  !identical(
    metadata$Source_Scrublet_Class,
    metadata$Source_Predicted_Doublet
  )) {
  stop("PRJCA015229 annotated human RNA/doublet design differs from audit")
}
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "PRJCA015229"
metadata$Technology <- ifelse(
  metadata$Tech == "scMultiome",
  "10x Genomics Multiome",
  "10x Genomics 3-prime Gene Expression"
)
metadata$Sequence <- "snRNA-seq"
metadata$Assay_Type <- ifelse(
  metadata$Tech == "scMultiome",
  "RNA arm of single-nucleus multiome",
  "single-nucleus RNA-seq"
)
metadata$Sequencing_Platform <- "Illumina NovaSeq 6000"
metadata$Library_Chemistry <- ifelse(
  metadata$Tech == "scMultiome",
  "10x Chromium Single Cell Multiome ATAC plus Gene Expression",
  "10x Chromium Single Cell 3-prime Gene Expression"
)
metadata$Source_RNA_Library_ID <- paste(
  metadata$Sample,
  metadata$Tech,
  sep = "::"
)
source_library_design <- unique(metadata[, c("Sample", "Tech")])
if (nrow(source_library_design) != 6L ||
  anyDuplicated(source_library_design$Sample) ||
  length(unique(metadata$Source_RNA_Library_ID)) != 6L) {
  stop("PRJCA015229 must retain one RNA library per human donor")
}
metadata$Sample_ID <- metadata$Source_RNA_Library_ID
metadata$Original_Donor_ID <- metadata$Sample
metadata$Original_Specimen_ID <- metadata$Sample
metadata$Original_Source_Record_ID <- metadata$Source_RNA_Library_ID
metadata$CellType_raw <- metadata$Annotation
metadata$Brain_Region <- "Anterior cingulate cortex"
metadata$Region <- "Anterior cingulate cortex"
metadata$Original_File_Type <- "10x Genomics filtered feature-barcode matrix"
metadata$Expression_Units <- "raw integer UMI counts"
metadata$Genome_Transcriptome_Build <- "GRCh38"
metadata$Raw_Integer_Counts_Available <- TRUE
metadata$Source_Full_Human_Annotated_RNA_Nuclei <- 40748L
metadata$Source_Retained_Human_Annotated_RNA_Nuclei <- 40748L
metadata$Source_Analysis_Singlet_Nuclei <- 40172L
metadata$Source_Preprocessing_Retained_Object_Scope <- paste(
  "complete publication-annotated human RNA nuclei from four Multiome",
  "RNA arms and two standalone snRNA-seq libraries; doublets retained"
)
metadata$Source_Preprocessing_Attrition_Status <- paste(
  "57,566 nuclei in the six human filtered feature-barcode matrices;",
  "40,748 have matching released publication annotations; all annotated",
  "nuclei are retained, including 576 Scrublet doublets that are excluded",
  "only from the primary analysis role"
)
metadata$Source_Preprocessing_Diagnosis_Audit_Status <- paste(
  "six neurologically normal adult male postmortem anterior cingulate",
  "cortex donors retained according to the source publication sample table"
)

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
sample_info$Sex_Source_Raw <- sample_info$Sex
sample_info$Sex_Source_Standardized <- sample_info$Sex
sample_info$Sex_Source_Column <- "source publication sample table"
sample_info$Sex_Assignment_Method <-
  "curated from source publication sample table; not inferred"

metadata <- merge(
  metadata, sample_info,
  by.x = "Sample", by.y = "Sample", all.x = TRUE
)
rownames(metadata) <- metadata$Cells

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "Assay_Type", "Sequencing_Platform",
  "Library_Chemistry", "Source_RNA_Library_ID",
  "Original_Donor_ID", "Original_Specimen_ID",
  "Original_Source_Record_ID", "CellType_raw", "BigCellType",
  "Annotation", "Source_Scrublet_Score", "Source_Scrublet_Class",
  "Source_Predicted_Doublet", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)

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
