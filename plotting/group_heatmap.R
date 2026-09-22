#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(scop)
  library(SeuratObject)
  library(ggplot2)
  library(ComplexHeatmap)
  library(grid)
})
source("functions/utils.R")
setDTthreads(2)
grDevices::pdfFonts(Arial = grDevices::pdfFonts("ArialMT")[[1]])
figure_dir <- "figures"
output_dir <- "results/annotation"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

formal <- as.data.table(read_celltype_assignments())
mapping <- unique(formal[, .(Cluster, CellType)])
type_counts <- formal[, .(Cells = .N), by = .(Group = CellType)]
stopifnot(
  nrow(mapping) == 75L,
  !anyDuplicated(mapping$Cluster),
  nrow(type_counts) == 16L,
  sum(type_counts$Cells) == 2712452L,
  type_counts[Group == "Excitatory neurons", Cells] == 971364L
)
type_levels <- names(brainomics_celltype_colors)[
  names(brainomics_celltype_colors) %in% type_counts$Group
]
rm(formal)
gc()

spec <- fread("results/annotation/display_markers.tsv")
stopifnot(
  nrow(spec) == 47L,
  !anyDuplicated(spec$Gene),
  setequal(spec$Marker_Group, type_levels)
)

evidence <- fread("results/annotation/full_cluster_marker_evidence.tsv")
evidence <- evidence[Gene %in% spec$Gene]
stopifnot(
  nrow(evidence) == 75L * nrow(spec),
  all(evidence$Available_Cells > 0),
  !anyDuplicated(evidence[, .(Cluster, Gene)])
)
evidence[, Group := mapping$CellType[match(Cluster, mapping$Cluster)]]
summary <- evidence[,
  .(
    Available_Cells = sum(Available_Cells),
    MeanLog = sum(MeanLog * Available_Cells) / sum(Available_Cells),
    PctPositive = sum(PctPositive * Available_Cells) / sum(Available_Cells)
  ),
  by = .(Group, Gene)
]
summary[, Cells := type_counts$Cells[match(Group, type_counts$Group)]]
stopifnot(
  nrow(summary) == 16L * nrow(spec),
  all(summary$Available_Cells <= summary$Cells)
)
fwrite(summary, file.path(output_dir, "celltype_marker_source.tsv"), sep = "\t")

grDevices::pdf(file = NULL, family = "Arial")
p <- celltype_summary_group_heatmap(
  as.data.frame(summary),
  as.data.frame(spec),
  type_levels
)
grDevices::dev.off()
ggsave(
  file.path(figure_dir, "group_heatmap.pdf"),
  p$plot,
  device = grDevices::cairo_pdf,
  width = 160,
  height = 160,
  units = "mm",
  family = "Arial",
  bg = "white"
)
message("Cell-type marker heatmap completed")

# Transposed companion reserves canvas width for full cell-type row labels.
grDevices::pdf(file = NULL, family = "Arial")
p_flip <- celltype_summary_group_heatmap(as.data.frame(summary), as.data.frame(spec),
  type_levels, flip = TRUE)
grDevices::dev.off()
ggsave(file.path(figure_dir, "group_heatmap_flip.pdf"), p_flip$plot,
  device = grDevices::cairo_pdf, width = 240, height = 160, units = "mm",
  family = "Arial", bg = "white")

# Cluster-level companion: retain every cluster instead of pooling cell types.
cluster_order <- mapping[order(match(CellType, type_levels), Cluster), Cluster]
cluster_type <- factor(mapping$CellType[match(cluster_order, mapping$Cluster)],
  levels = type_levels)
genes <- spec$Gene
idx <- match(paste(rep(cluster_order, each = length(genes)),
  rep(genes, length(cluster_order))), paste(evidence$Cluster, evidence$Gene))
stopifnot(!anyNA(idx), length(unique(cluster_order)) == 75L)
cluster_means <- matrix(evidence$MeanLog[idx], nrow = length(genes),
  dimnames = list(genes, cluster_order))
cluster_fraction <- matrix(evidence$PctPositive[idx], nrow = length(genes),
  dimnames = dimnames(cluster_means))
cluster_z <- t(scale(t(cluster_means)))
stopifnot(all(is.finite(cluster_z)), all(is.finite(cluster_fraction)),
  all(cluster_fraction >= 0 & cluster_fraction <= 1))
marker_levels <- unique(spec$Marker_Group)
cluster_col <- circlize::colorRamp2(c(-2, -1, 0, 1, 2),
  c("#3288BD", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D53E4F"))
legend_text <- gpar(fontfamily = "Arial", fontsize = 7)
cluster_heatmap <- Heatmap(cluster_z,
  name = "Z-score",
  col = cluster_col, cluster_rows = FALSE, cluster_columns = FALSE,
  cluster_row_slices = FALSE, cluster_column_slices = FALSE,
  row_split = factor(spec$Marker_Group, levels = marker_levels),
  column_split = cluster_type, row_title = NULL, column_title = NULL,
  row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
  row_names_side = "right", row_names_gp = gpar(fontfamily = "Arial",
    fontsize = 7, fontface = "italic"),
  column_names_rot = 90, column_names_gp = legend_text,
  rect_gp = gpar(type = "none"), border = "grey75", use_raster = FALSE,
  width = unit(220, "mm"), height = unit(158, "mm"),
  top_annotation = HeatmapAnnotation(Celltype = cluster_type,
    col = list(Celltype = brainomics_celltype_colors[type_levels]),
    show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"),
    annotation_legend_param = list(Celltype = list(
      title_gp = legend_text, labels_gp = legend_text))),
  left_annotation = rowAnnotation(Marker = factor(spec$Marker_Group,
      levels = marker_levels),
    col = list(Marker = brainomics_celltype_colors[marker_levels]),
    show_annotation_name = FALSE, show_legend = FALSE,
    simple_anno_size = unit(2, "mm")),
  heatmap_legend_param = list(title = "Z-score\n(mean log expression)",
    at = c(-2, -1, 0, 1, 2), labels = c("\u2264 -2", "-1", "0", "1", "\u2265 2"),
    title_gp = legend_text, labels_gp = legend_text),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (cluster_fraction[i, j] > 0) grid.points(x, y, pch = 21,
      size = unit(2.2 * cluster_fraction[i, j], "mm"),
      gp = gpar(fill = fill, col = "grey20", lwd = .35))
  })
cluster_percent_legend <- Legend(title = "Percent detected",
  labels = paste0(seq(20, 100, 20), "%"), type = "points", pch = 21,
  size = unit(2.2 * seq(.2, 1, .2), "mm"),
  legend_gp = gpar(fill = "grey60", col = "grey20", lwd = .35),
  title_gp = legend_text, labels_gp = legend_text)
grDevices::cairo_pdf(file.path(figure_dir, "group_heatmap_clusters.pdf"),
  width = 300 / 25.4, height = 185 / 25.4, family = "Arial")
draw(cluster_heatmap, heatmap_legend_side = "right",
  annotation_legend_side = "right", merge_legends = TRUE,
  heatmap_legend_list = list(cluster_percent_legend),
  padding = unit(c(4, 4, 4, 4), "mm"))
grDevices::dev.off()
message("75-cluster marker heatmap completed")

# Optional layout trial: selected identity contrasts, not measured validation calls.
# Human brain marker basis: https://doi.org/10.1038/s41591-024-03150-z
# The +/- contrasts below are interpretation aids for canonical neuronal identity;
# they are not universal absence rules or MGE/CGE lineage validation.
if ("--marker-relations-preview" %in% commandArgs(trailingOnly = TRUE)) {
  preview <- cluster_heatmap
  preview@matrix_param$row_gap <- unit(3, "mm")
  grDevices::cairo_pdf(file.path(figure_dir,
    "group_heatmap_clusters_marker_relations_preview.pdf"),
    width = 300 / 25.4, height = 200 / 25.4, family = "Arial")
  draw(preview, heatmap_legend_side = "right",
    annotation_legend_side = "right", merge_legends = TRUE,
    heatmap_legend_list = list(cluster_percent_legend),
    padding = unit(c(13, 4, 4, 4), "mm"))
  contrasts <- data.frame(
    MarkerGroup = rep(c("Excitatory neurons", "MGE-derived inhibitory neurons"), each = 3),
    Celltype = rep(c("Excitatory neurons", "MGE-derived inhibitory neurons",
      "CGE-derived inhibitory neurons"), 2),
    Sign = c("+", "\u2212", "\u2212", "\u2212", "+", "+"))
  for (k in seq_len(nrow(contrasts))) {
    r <- match(contrasts$MarkerGroup[k], marker_levels)
    c <- match(contrasts$Celltype[k], type_levels)
    # SLC17A7 is row 3/3; GAD1/GAD2 are rows 1/3 and 2/3 (exclude NXPH1).
    glutamate <- contrasts$MarkerGroup[k] == "Excitatory neurons"
    decorate_heatmap_body("Z-score", {
      grid.rect(x = .5, y = if (glutamate) 1/6 else 2/3,
        width = unit(1, "npc") - unit(.3, "mm"),
        height = unit(if (glutamate) 1/3 else 2/3, "npc"),
        gp = gpar(fill = NA, col = "grey15", lwd = .7, lty = 2))
      grid.text(contrasts$Sign[k], x = if (glutamate) .25 else .75,
        y = if (glutamate) unit(-1.3, "mm") else unit(1, "npc") + unit(1.3, "mm"),
        gp = gpar(fontfamily = "Arial", fontsize = 8))
    }, row_slice = r, column_slice = c)
  }
  grid.text(paste0("Selected neuronal identity contrasts:  + supporting marker;  \u2212 competing-identity marker.\n",
    "Symbols indicate expected relationships, not measured positivity/negativity; unmarked combinations are not assessed."),
    x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"),
    gp = gpar(fontfamily = "Arial", fontsize = 7))
  grDevices::dev.off()
  message("Marker-relation layout preview completed")
}
