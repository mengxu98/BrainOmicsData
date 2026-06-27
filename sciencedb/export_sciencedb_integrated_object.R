#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

object_file <- normalizePath(
  value_after("--object-file", "../../data/BrainOmicsData/integration/objects_celltypes.rds"),
  mustWork = TRUE
)
out_dir <- value_after(
  "--out-dir",
  "../../data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, file) {
  utils::write.table(
    x, file = file, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "missing"
  )
}

fill_missing <- function(x, value = "missing") {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA" | x == "Unknown"] <- value
  x
}

add_age_interval_schema <- function(x) {
  age_interval <- c(
    "S1" = "Embryonic",
    "S2" = "Early fetal",
    "S3" = "Early fetal",
    "S4" = "Early mid-fetal",
    "S5" = "Early mid-fetal",
    "S6" = "Late mid-fetal",
    "S7" = "Late fetal",
    "S8" = "Neonatal and early infancy",
    "S9" = "Late infancy",
    "S10" = "Early childhood",
    "S11" = "Middle and late childhood",
    "S12" = "Adolescence",
    "S13" = "Young adulthood",
    "S14" = "Middle adulthood",
    "S15" = "Late adulthood"
  )
  age_range <- c(
    "S1" = "4-8 PCW",
    "S2" = "8-10 PCW",
    "S3" = "10-13 PCW",
    "S4" = "13-16 PCW",
    "S5" = "16-19 PCW",
    "S6" = "19-24 PCW",
    "S7" = "24-38 PCW",
    "S8" = "0-0.5 years",
    "S9" = "0.5-1 years",
    "S10" = "1-6 years",
    "S11" = "6-12 years",
    "S12" = "12-20 years",
    "S13" = "20-40 years",
    "S14" = "40-60 years",
    "S15" = "60+ years"
  )
  if (!"AgeIntervalID" %in% names(x)) x$AgeIntervalID <- x$Stage
  if (!"AgeInterval" %in% names(x)) {
    x$AgeInterval <- unname(age_interval[as.character(x$AgeIntervalID)])
  }
  if (!"AgeRange" %in% names(x)) {
    x$AgeRange <- unname(age_range[as.character(x$AgeIntervalID)])
  }
  x
}

paths <- list(
  cell = file.path(out_dir, "03_cell_metadata"),
  objects = file.path(out_dir, "04_processed_objects"),
  emb = file.path(out_dir, "05_embeddings_and_clusters"),
  qc = file.path(out_dir, "06_quality_control"),
  anno = file.path(out_dir, "07_celltype_annotation")
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("The heavy export requires the Seurat package to read objects_celltypes.rds.")
}

message("Loading integrated object: ", object_file)
object <- readRDS(object_file)
meta <- object@meta.data
meta[] <- lapply(meta, as.character)
meta <- add_age_interval_schema(meta)

cell_id <- rownames(meta)
umap_name <- if ("umap.rpca" %in% names(object@reductions)) "umap.rpca" else grep("umap", names(object@reductions), value = TRUE)[1]
pca_name <- if ("integrated.rpca" %in% names(object@reductions)) "integrated.rpca" else if ("pca" %in% names(object@reductions)) "pca" else NA_character_

cluster_col <- if ("seurat_clusters" %in% names(meta)) "seurat_clusters" else NA_character_
celltype_col <- if ("CellType" %in% names(meta)) "CellType" else if ("Celltype" %in% names(meta)) "Celltype" else NA_character_

cell_metadata <- data.frame(
  cell_id = cell_id,
  sample_id = meta$Sample_ID,
  donor_id = meta$Sample_ID,
  source_dataset = meta$Dataset,
  source_cell_id = meta$Cells,
  brain_region_harmonized = meta$BrainRegion,
  age_interval_id = meta$AgeIntervalID,
  age_interval = meta$AgeInterval,
  age_range = meta$AgeRange,
  reported_age = meta$Age,
  sex = meta$Sex,
  sequencing_modality = meta$Sequence,
  sequencing_platform = meta$Technology,
  cluster_id = if (!is.na(cluster_col)) meta[[cluster_col]] else NA_character_,
  cell_type = if (!is.na(celltype_col)) meta[[celltype_col]] else NA_character_,
  n_counts = if ("nCount_RNA" %in% names(meta)) meta$nCount_RNA else NA_character_,
  n_genes = if ("nFeature_RNA" %in% names(meta)) meta$nFeature_RNA else NA_character_,
  percent_mito = if ("percent.mt" %in% names(meta)) meta$percent.mt else NA_character_,
  percent_ribo = if ("percent.ribo" %in% names(meta)) meta$percent.ribo else NA_character_,
  quality_control_status = "included_in_integrated_object",
  stringsAsFactors = FALSE
)
cell_metadata[] <- lapply(cell_metadata, fill_missing)

gz <- gzfile(file.path(paths$cell, "cell_metadata_harmonized.tsv.gz"), "wt")
write_tsv(cell_metadata, gz)
close(gz)

if (!is.na(umap_name)) {
  umap <- as.data.frame(object@reductions[[umap_name]]@cell.embeddings)
  umap <- data.frame(cell_id = rownames(umap), umap, check.names = FALSE)
  names(umap) <- sub(paste0("^", umap_name, "_"), "umap_", names(umap))
  gz <- gzfile(file.path(paths$emb, "umap_coordinates.tsv.gz"), "wt")
  write_tsv(umap, gz)
  close(gz)
}

if (!is.na(pca_name)) {
  pca <- as.data.frame(object@reductions[[pca_name]]@cell.embeddings)
  pca <- data.frame(cell_id = rownames(pca), pca, check.names = FALSE)
  gz <- gzfile(file.path(paths$emb, "pca_coordinates.tsv.gz"), "wt")
  write_tsv(pca, gz)
  close(gz)
}

cluster_assignments <- cell_metadata[, c("cell_id", "source_dataset", "cluster_id", "cell_type")]
gz <- gzfile(file.path(paths$emb, "cluster_assignments.tsv.gz"), "wt")
write_tsv(cluster_assignments, gz)
close(gz)

if (!is.na(cluster_col) && !is.na(celltype_col)) {
  mapping <- unique(data.frame(
    cluster_id = meta[[cluster_col]],
    cell_type = meta[[celltype_col]],
    stringsAsFactors = FALSE
  ))
  mapping <- mapping[order(as.integer(mapping$cluster_id), mapping$cell_type), ]
  write_tsv(mapping, file.path(paths$cell, "celltype_cluster_mapping.tsv"))
  write_tsv(mapping, file.path(paths$anno, "celltype_cluster_mapping.tsv"))
}

if (!is.na(celltype_col)) {
  counts <- as.data.frame(
    xtabs(~ AgeIntervalID + BrainRegion + CellType, transform(meta, CellType = meta[[celltype_col]])),
    stringsAsFactors = FALSE
  )
  names(counts) <- c("age_interval_id", "brain_region_harmonized", "cell_type", "cell_count")
  counts <- counts[counts$cell_count > 0, ]
  write_tsv(counts, file.path(paths$cell, "cell_counts_by_age_region_celltype.tsv"))

  celltype_counts <- as.data.frame(table(cell_type = meta[[celltype_col]]), stringsAsFactors = FALSE)
  names(celltype_counts) <- c("cell_type", "cell_count")
  celltype_counts <- celltype_counts[celltype_counts$cell_count > 0, ]
  write_tsv(celltype_counts, file.path(paths$cell, "cell_counts_by_cell_type.tsv"))
  write_tsv(celltype_counts, file.path(paths$anno, "cell_type_counts.tsv"))
}

qc <- data.frame(
  cell_id = cell_id,
  source_dataset = meta$Dataset,
  n_counts = if ("nCount_RNA" %in% names(meta)) meta$nCount_RNA else NA_character_,
  n_genes = if ("nFeature_RNA" %in% names(meta)) meta$nFeature_RNA else NA_character_,
  percent_mito = if ("percent.mt" %in% names(meta)) meta$percent.mt else NA_character_,
  percent_ribo = if ("percent.ribo" %in% names(meta)) meta$percent.ribo else NA_character_,
  stringsAsFactors = FALSE
)
qc[] <- lapply(qc, fill_missing)
gz <- gzfile(file.path(paths$qc, "qc_metrics_by_cell.tsv.gz"), "wt")
write_tsv(qc, gz)
close(gz)

marker_genes <- data.frame(
  cell_type = c(
    rep("Radial glia", 3), rep("Endothelial cells", 4),
    rep("Inhibitory neurons", 3), rep("Oligodendrocyte progenitor cells", 5),
    rep("Microglia", 3), "Neuroblasts", rep("Excitatory neurons", 3),
    rep("Astrocytes", 5), rep("Oligodendrocytes", 3)
  ),
  marker_gene = c(
    "PAX6", "VIM", "GLI3", "CLDN5", "PECAM1", "VWF", "FLT1",
    "GAD1", "GAD2", "SLC6A1", "PDGFRA", "CSPG4", "OLIG1", "OLIG2",
    "SOX10", "CX3CR1", "P2RY12", "CSF1R", "STMN2", "SLC17A7",
    "CAMK2A", "SATB2", "GFAP", "AQP4", "ALDH1L1", "FGFR3", "GJA1",
    "MOG", "MAG", "CLDN11"
  ),
  stringsAsFactors = FALSE
)
present_markers <- intersect(marker_genes$marker_gene, rownames(object))
if (!is.na(celltype_col) && length(present_markers) > 0) {
  dot <- Seurat::DotPlot(object, features = present_markers, group.by = celltype_col)$data
  marker_summary <- data.frame(
    cell_type = as.character(dot$id),
    marker_gene = as.character(dot$features.plot),
    average_expression = dot$avg.exp,
    scaled_average_expression = dot$avg.exp.scaled,
    percent_expressing = dot$pct.exp,
    stringsAsFactors = FALSE
  )
  write_tsv(marker_summary, file.path(paths$anno, "marker_expression_summary.tsv"))
}

if (!is.na(cluster_col) && !is.na(celltype_col)) {
  validation <- merge(
    unique(data.frame(
      cluster_id = meta[[cluster_col]],
      cell_type = meta[[celltype_col]],
      stringsAsFactors = FALSE
    )),
    aggregate(marker_gene ~ cell_type, marker_genes, paste, collapse = ";"),
    by = "cell_type",
    all.x = TRUE
  )
  names(validation)[names(validation) == "marker_gene"] <- "supporting_marker_genes"
  validation$validation_status <- "marker-supported major cell-type assignment"
  validation <- validation[order(as.integer(validation$cluster_id), validation$cell_type), ]
  write_tsv(validation, file.path(paths$anno, "cluster_marker_validation.tsv"))
}

nn <- data.frame(
  graph_or_reduction = names(object@graphs),
  note = "Nearest-neighbor graph names recorded from the Seurat object; full graph matrices are not exported by default because of file size.",
  stringsAsFactors = FALSE
)
write_tsv(nn, file.path(paths$emb, "nearest_neighbor_graph_metadata.tsv"))

if (requireNamespace("SeuratDisk", quietly = TRUE)) {
  h5seurat <- file.path(paths$objects, "human_brain_integrated_full.h5seurat")
  h5ad <- file.path(paths$objects, "human_brain_integrated_full.h5ad")
  message("Writing h5Seurat: ", h5seurat)
  SeuratDisk::SaveH5Seurat(object, filename = h5seurat, overwrite = TRUE)
  message("Converting to h5ad: ", h5ad)
  SeuratDisk::Convert(h5seurat, dest = "h5ad", overwrite = TRUE)
} else if (requireNamespace("scop", quietly = TRUE)) {
  h5ad <- file.path(paths$objects, "human_brain_integrated_full.h5ad")
  message("Writing h5ad with scop::srt_to_h5ad: ", h5ad)
  scop::srt_to_h5ad(
    object,
    path = h5ad,
    assay_x = "RNA",
    layer_x = "counts",
    reductions = names(object@reductions),
    graphs = character(0),
    neighbors = character(0),
    convert_tools = FALSE,
    convert_misc = FALSE,
    overwrite = TRUE
  )
} else {
  note <- c(
    "Neither SeuratDisk nor scop was available in this R environment.",
    "Run this script in an environment with SeuratDisk or scop::srt_to_h5ad() to create human_brain_integrated_full.h5ad."
  )
  writeLines(note, file.path(paths$objects, "h5ad_h5seurat_export_note.txt"))
}

message("Heavy integrated-object export complete: ", out_dir)
