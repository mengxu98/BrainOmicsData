#!/usr/bin/env Rscript

# RPCA cluster and cell-type UMAPs for the annotation row.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(scop)
  library(Seurat)
})

source("functions/utils.R")

source("functions/config.R")
panel_dir <- "figures";dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE);four_method_only<-FALSE
metadata <- fread(resolve_input_file(
  "BRAINOMICS_FIG2_METADATA_FILE", root, "core_metadata_minimal.tsv.gz"
))
embedding <- readRDS(resolve_input_file(
  "BRAINOMICS_RPCA_UMAP_FILE", root, "embedding_umap.rpca.rds"
))
stopifnot(is.matrix(embedding),nrow(embedding)==2602031L,identical(rownames(embedding),metadata$Cells))
umap <- data.frame(Cell=metadata$Cells,RPCA_1=embedding[,1],RPCA_2=embedding[,2])
assignments <- data.frame(Cells=metadata$Cells,Cluster=metadata$Cluster,CellType=a$Working_CellType[match(metadata$Cluster,a$Cluster)])
stopifnot(!anyNA(assignments$CellType),length(unique(assignments$Cluster))==75L)
display_celltype <- assignments$CellType
celltype_order <- names(brainomics_celltype_colors)
celltype_colors <- brainomics_celltype_colors

cells <- as.character(umap$Cell)
# Sort cluster labels numerically while preserving their existing colors.
cluster_color_order <- sort(unique(as.character(assignments$Cluster)))
stopifnot(all(grepl("^C[0-9]+$", cluster_color_order)))
cluster_order <- cluster_color_order[order(as.integer(sub("^C", "", cluster_color_order)))]
counts <- sparseMatrix(
  i = integer(), j = integer(), x = numeric(),
  dims = c(2L, length(cells)),
  dimnames = list(c("placeholderA", "placeholderB"), cells)
)
plot_metadata <- data.frame(
  CellType = factor(display_celltype, levels = celltype_order),
  Cluster = factor(assignments$Cluster, levels = cluster_order),
  row.names = cells
)
object <- CreateSeuratObject(
  counts = methods::as(counts, "dgCMatrix"), meta.data = plot_metadata,
  min.cells = 0L, min.features = 0L
)
rpca_embedding <- as.matrix(umap[, c("RPCA_1", "RPCA_2")])
rownames(rpca_embedding) <- cells
colnames(rpca_embedding) <- c("UMAP_1", "UMAP_2")
object[["umap.rpca"]] <- CreateDimReducObject(
  embeddings = rpca_embedding, key = "UMAP_", assay = "RNA"
)

set.seed(20260906)
draw_order <- colnames(object)[sample.int(ncol(object))]
if (!four_method_only) {
  p_celltype <- scop::CellDimPlot(
    object,
    reduction = "umap.rpca", group.by = "CellType",
    palcolor = celltype_colors, label = FALSE, seed = 11,
    raster = TRUE, raster.dpi = c(1600, 1600), pt.size = 2,
    xlab = "UMAP_1", ylab = "UMAP_2", legend.title = "Cell type",
    theme_use = "theme_blank_axis", theme_args = list(text = element_text(family = "Arial")),
    combine = FALSE
  )[[1]]
  p_celltype <- order_dim_plot_cells(p_celltype, draw_order)
  stopifnot(identical(rownames(p_celltype$data), draw_order))
  message("RPCA cell-type UMAP completed: ", panel_dir)
}

if (!four_method_only) {
  cluster_levels <- levels(object$Cluster)
  cluster_colors <- setNames(grDevices::hcl.colors(length(cluster_color_order), "Dynamic"), cluster_color_order)
  cluster_colors <- cluster_colors[cluster_levels]
  p_cluster <- scop::CellDimPlot(object, reduction = "umap.rpca", group.by = "Cluster",
    palcolor = cluster_colors, label = FALSE, seed = 11, show_stat = FALSE,
    raster = TRUE, raster.dpi = c(1600, 1600), pt.size = 2,
    xlab = "UMAP_1", ylab = "UMAP_2", legend.title = "Cluster",
    theme_use = "theme_blank_axis", theme_args = list(text = element_text(family = "Arial")),
    combine = FALSE)[[1]]
  p_cluster <- order_dim_plot_cells(p_cluster, draw_order)
  # Keep both legends beside their UMAPs at manuscript width.
  annotation_theme <- theme(
    text = element_text(family = "Arial", size = 6, face = "plain"),
    plot.title = element_text(size = 6.5, face = "plain"),
    plot.subtitle = element_text(size = 6), axis.title = element_blank(),
    legend.position = "right", legend.text = element_text(size = 5),
    legend.title = element_text(size = 5.5),
    legend.key.height = grid::unit(2, "mm"), legend.key.width = grid::unit(2, "mm"),
    legend.spacing.x = grid::unit(.5, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "mm"))
  p_cluster <- p_cluster + labs(title = "RPCA clusters (n = 75)", subtitle = NULL) +
    annotation_theme + guides(colour = guide_legend(ncol = 3, byrow = FALSE,
      override.aes = list(size = 1.1)))
  p_annotated <- p_celltype + labs(title = "Cell types (n = 12)", subtitle = NULL) +
    annotation_theme + guides(colour = guide_legend(ncol = 1, byrow = TRUE,
      override.aes = list(size = 1.1)))
  # Reapply the shared arrow axes at the final panel font size.
  annotation_plots <- lapply(list(p_cluster, p_annotated), function(p) {
    p$layers <- Filter(function(layer) !inherits(layer$geom, "GeomCustomAnn"), p$layers)
    p + theme_blank_axis(lab_size = 5.5, axis_lwd = .6) + p$theme +
      p$coordinates + annotation_theme
  })
  annotation_row <- patchwork::wrap_plots(annotation_plots, nrow = 1)
  ggplot2::ggsave(file.path(panel_dir, "fig2f.pdf"), annotation_row,
    device = grDevices::cairo_pdf, width = 168, height = 66, units = "mm",
    family = "Arial", bg = "white")
}
