source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/integration/"
fig_dir <- check_dir("figures/")

objects_plot <- readRDS(
  file.path(data_dir, "objects_celltype_plot.rds")
)


marker_genes <- c(
  # Radial glia
  "PAX6", "VIM", "GLI3",
  # Endothelial cells
  "CLDN5", "PECAM1", "VWF", "FLT1",
  # Inhibitory neurons
  "GAD1", "GAD2", "SLC6A1",
  # Oligodendrocyte progenitor cells (OPCs)
  "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10",
  # Microglia
  "CX3CR1", "P2RY12", "CSF1R",
  # Neuroblasts
  "STMN2",
  # Excitatory neurons
  "SLC17A7", "CAMK2A", "SATB2",
  # Astrocytes
  "GFAP", "AQP4", "ALDH1L1", "FGFR3", "GJA1",
  # Oligodendrocytes
  "MOG", "MAG", "CLDN11"
)

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

objects_plot$Stage <- factor(
  objects_plot$Stage,
  levels = paste0("S", 1:15)
)
objects_plot$AgeIntervalID <- factor(
  objects_plot$AgeIntervalID,
  levels = paste0("S", 1:15)
)

Idents(objects_plot) <- "CellType"
celltype_levels <- sort(as.character(unique(objects_plot$CellType)))
objects_plot$Celltype <- factor(
  objects_plot$CellType,
  levels = celltype_levels
)

development_stages <- c(
  "Embryonic",
  "Early fetal",
  "Early fetal",
  "Early mid-fetal",
  "Early mid-fetal",
  "Late mid-fetal",
  "Late fetal",
  "Neonatal and early infancy",
  "Late infancy",
  "Early childhood",
  "Middle and late childhood",
  "Adolescence",
  "Young adulthood",
  "Middle adulthood",
  "Late adulthood"
)
development_stages2 <- c(
  "Embryonic (4-8 PCW)",
  "Early fetal (8-10 PCW)",
  "Early fetal (10-13 PCW)",
  "Early mid-fetal (13-16 PCW)",
  "Early mid-fetal (16-19 PCW)",
  "Late mid-fetal (19-24 PCW)",
  "Late fetal (24-38 PCW)",
  "Neonatal and early infancy (0-0.5Y)",
  "Late infancy (0.5-1Y)",
  "Early childhood (1-6Y)",
  "Middle and late childhood (6-12Y)",
  "Adolescence (12-20Y)",
  "Young adulthood (20-40Y)",
  "Middle adulthood (40-60Y)",
  "Late adulthood (60+Y)"
)
objects_plot$DevelopmentStage <- factor(
  objects_plot$DevelopmentStage,
  levels = development_stages2
)
objects_plot$AgeInterval <- factor(
  objects_plot$AgeInterval,
  levels = unique(development_stages)
)

stage_def <- data.frame(
  AgeIntervalID = paste0("S", 1:15),
  AgeInterval = development_stages,
  AgeRange = c(
    "4-8 PCW", "8-10 PCW", "10-13 PCW", "13-16 PCW", "16-19 PCW",
    "19-24 PCW", "24-38 PCW",
    "0-0.5 years", "0.5-1 years", "1-6 years", "6-12 years",
    "12-20 years", "20-40 years", "40-60 years", "60+ years"
  ),
  stringsAsFactors = FALSE
)
stage_counts <- as.data.frame(
  table(AgeIntervalID = objects_plot$AgeIntervalID),
  stringsAsFactors = FALSE
)

names(stage_counts)[2] <- "N_cells"
stage_def <- merge(stage_def, stage_counts, by = "AgeIntervalID", all.x = TRUE)
stage_def$N_cells[is.na(stage_def$N_cells)] <- 0
stage_def$AgeIntervalID <- factor(stage_def$AgeIntervalID, levels = paste0("S", 1:15))
stage_def <- stage_def[order(stage_def$AgeIntervalID), ]
stage_def$AgeIntervalID <- as.character(stage_def$AgeIntervalID)

stage_labels <- as.matrix(
  stage_def[, c("AgeIntervalID", "AgeInterval", "AgeRange", "N_cells")]
)
stage_labels[, "N_cells"] <- format(
  stage_def$N_cells,
  big.mark = ",", trim = TRUE
)
mat_text <- matrix(0, nrow = 15, ncol = 4)
colnames(mat_text) <- c(
  "Age interval ID", "Age interval", "Age range", "Cell count"
)
rownames(mat_text) <- as.character(stage_def$AgeIntervalID)

stage_table_gp <- gpar(fontsize = 9)
stage_table_padding_mm <- 1.5
stage_cell_fill <- "#EFF7FC"
stage_cell_border <- "#9fbbd3ff"
stage_text_width_mm <- function(labels) {
  label_widths <- lapply(
    labels,
    function(label) grobWidth(textGrob(label, gp = stage_table_gp))
  )
  convertWidth(max(do.call(unit.c, label_widths)), "mm", TRUE)
}
stage_column_widths_mm <- vapply(
  seq_len(ncol(stage_labels)),
  function(j) {
    stage_text_width_mm(c(colnames(mat_text)[j], stage_labels[, j])) +
      stage_table_padding_mm * 2
  },
  numeric(1)
)
stage_color_width_mm <- 5
stage_row_height_mm <- 5.5
stage_header_height_mm <- 5
stage_table_width_mm <- sum(stage_column_widths_mm)
stage_plot_width_mm <- stage_color_width_mm + stage_table_width_mm
stage_plot_height_mm <- stage_row_height_mm * nrow(stage_labels) +
  stage_header_height_mm
stage_pdf_width <- (stage_plot_width_mm + 4) / 25.4
stage_pdf_height <- (stage_plot_height_mm + 2) / 25.4

draw_stage_table <- function() {
  grid.newpage()
  pushViewport(
    viewport(
      x = unit(2, "mm"), y = unit(1, "mm"),
      width = unit(stage_plot_width_mm, "mm"),
      height = unit(stage_plot_height_mm, "mm"),
      just = c("left", "bottom")
    )
  )

  n_rows <- nrow(stage_labels)
  table_left_mm <- stage_color_width_mm
  table_bottom_mm <- stage_header_height_mm

  grid.rect(
    x = unit(0, "mm"),
    y = unit(
      table_bottom_mm + stage_row_height_mm * (n_rows - seq_len(n_rows)),
      "mm"
    ),
    width = unit(stage_color_width_mm, "mm"),
    height = unit(stage_row_height_mm, "mm"),
    gp = gpar(fill = color_stages[as.character(stage_def$AgeIntervalID)], col = NA),
    just = c("left", "bottom")
  )

  col_left_mm <- table_left_mm
  for (j in seq_len(ncol(stage_labels))) {
    for (i in seq_len(n_rows)) {
      grid.rect(
        x = unit(col_left_mm, "mm"),
        y = unit(
          table_bottom_mm + stage_row_height_mm * (n_rows - i),
          "mm"
        ),
        width = unit(stage_column_widths_mm[j], "mm"),
        height = unit(stage_row_height_mm, "mm"),
        gp = gpar(fill = stage_cell_fill, col = stage_cell_border, lwd = 0.5),
        just = c("left", "bottom")
      )
      grid.text(
        stage_labels[i, j],
        x = unit(col_left_mm + stage_table_padding_mm, "mm"),
        y = unit(
          table_bottom_mm + stage_row_height_mm * (n_rows - i + 0.5),
          "mm"
        ),
        gp = stage_table_gp,
        just = "left"
      )
    }
    grid.text(
      colnames(mat_text)[j],
      x = unit(col_left_mm + stage_column_widths_mm[j] / 2, "mm"),
      y = unit(stage_header_height_mm / 2, "mm"),
      gp = stage_table_gp
    )
    col_left_mm <- col_left_mm + stage_column_widths_mm[j]
  }

  popViewport()
}

pdf(
  file.path(fig_dir, "development_stage_annotation.pdf"),
  width = stage_pdf_width, height = stage_pdf_height
)
draw_stage_table()
dev.off()

p1 <- DotPlot(
  objects_plot,
  features = marker_genes,
  group.by = "seurat_clusters",
  cols = c("gray80", "#15559A"),
  dot.scale = 5
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.text = element_text(size = 12)
  ) +
  coord_fixed()

ggsave(
  file.path(fig_dir, "dot_plot_markergenes.pdf"),
  p1,
  width = 9,
  height = 28
)

p2 <- FeatureDimPlot(
  objects_plot,
  features = marker_genes,
  reduction = "umap.rpca",
  ncol = 6,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)
ggsave(
  file.path(fig_dir, "feature_plots_rpca.pdf"),
  p2,
  width = 13,
  height = 10
)

p3 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "seurat_clusters",
  palcolor = colors128,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)
ggsave(
  file.path(fig_dir, "dim_plots_seurat_clusters_rpca.pdf"),
  p3,
  width = 11, height = 4.5
)

p4 <- CellDimPlot(
  objects_plot,
  reduction = "umap.unintegrated",
  group.by = "Dataset",
  palcolor = colors32,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)

p5 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "BrainRegion",
  palcolor = colors128,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)

color_stages2 <- color_stages[!duplicated(development_stages)]
names(color_stages2) <- unique(development_stages)
p6 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "AgeInterval",
  palcolor = color_stages2,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)
p7 <- p4 + p5 + p6
ggsave(
  file.path(fig_dir, "brainregion_stage_rpca.pdf"),
  p11,
  width = 28, height = 4
)

p8 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "Dataset",
  palcolor = colors32,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)
p9 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "AgeIntervalID",
  palcolor = color_stages,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)
p10 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "Celltype",
  palcolor = color_celltypes[celltype_levels],
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)

p11 <- p8 + p9 + p10
ggsave(
  file.path(fig_dir, "dim_plots_datasets_stage_celltype.pdf"),
  p13,
  width = 19, height = 4
)

marker_genes_split <- data.frame(
  Celltype = c(
    rep("Radial glia", 3),
    rep("Endothelial cells", 4),
    rep("Inhibitory neurons", 3),
    rep("Oligodendrocyte progenitor cells", 5),
    rep("Microglia", 3),
    rep("Neuroblasts", 1),
    rep("Excitatory neurons", 3),
    rep("Astrocytes", 5),
    rep("Oligodendrocytes", 3)
  ),
  Genes = marker_genes
)
celltype_levels <- sort(as.character(unique(marker_genes_split$Celltype)))
marker_genes_split$Celltype <- factor(
  marker_genes_split$Celltype,
  levels = celltype_levels
)
objects_plot$Celltype <- factor(
  objects_plot$Celltype,
  levels = celltype_levels
)
celltype_colors <- color_celltypes[celltype_levels]
gh <- GroupHeatmap(
  objects_plot,
  exp_legend_title = "Z-score",
  features = marker_genes_split$Genes,
  feature_split = marker_genes_split$Celltype,
  group.by = "Celltype",
  group_palcolor = celltype_colors,
  cell_annotation_palcolor = celltype_colors,
  feature_split_palcolor = celltype_colors,
  heatmap_palette = "Spectral",
  height = 6,
  width = 3,
  add_dot = TRUE,
  dot_size = unit(6, "mm"),
  nlabel = 0,
  show_row_names = TRUE,
  border = TRUE,
  ht_params = list(
    row_names_gp = gpar(fontface = "italic")
  )
)
pdf(
  file.path(fig_dir, "group_heatmap_markergenes.pdf"),
  width = 9.2,
  height = 6.7
)
print(gh$plot)
dev.off()


p_datasets <- p4 + p8 + p3
ggsave(
  file.path(fig_dir, "dim_plots_datasets_clusters.pdf"),
  p_datasets,
  width = 19, height = 4
)

p2_gh <- p2 / gh$plot

pdf(
  file.path(fig_dir, "dims_group_heatmap_markergenes.pdf"),
  width = 9.2,
  height = 15
)
print(p2_gh)
dev.off()
