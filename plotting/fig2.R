source("functions/prepare_env.R")
source("functions/metadata_schema.R")

data_dir <- "../../data/BrainOmicsData/integration/"
fig_dir <- check_dir("figures/")
res_dir <- check_dir("results/")

objects_plot <- readRDS(
  file.path(data_dir, "objects_celltype_plot.rds")
)
objects_plot@meta.data <- add_metadata_schema(objects_plot@meta.data)

if (!file.exists(file.path(res_dir, "lisi_results.rds"))) {
  lisi_data <- readRDS(file.path(data_dir, "lisi_data.rds"))
  raw <- compute_lisi(
    X = lisi_data[["umap_raw"]],
    meta_data = lisi_data[["meta_data"]],
    label_colnames = "Dataset"
  )
  harmony_lisi <- compute_lisi(
    X = lisi_data[["umap_harmony"]],
    meta_data = lisi_data[["meta_data"]],
    label_colnames = "Dataset"
  )
  rpca_lisi <- compute_lisi(
    X = lisi_data[["umap_rpca"]],
    meta_data = lisi_data[["meta_data"]],
    label_colnames = "Dataset"
  )
  lisi_results <- cbind(
    raw,
    harmony_lisi,
    rpca_lisi
  )
  names(lisi_results) <- c("Raw", "RPCA")

  saveRDS(lisi_results, file.path(res_dir, "lisi_results.rds"))
} else {
  lisi_results <- readRDS(file.path(res_dir, "lisi_results.rds"))
}

lisi_long <- tidyr::gather(
  lisi_results,
  key = "Method",
  value = "LISI"
)
lisi_long$Method <- factor(
  lisi_long$Method,
  levels = c("Raw", "RPCA")
)

lisi_long_2 <- lisi_long[lisi_long$Method %in% c("Raw", "RPCA"), ]
lisi_long_2$Method <- factor(
  lisi_long_2$Method,
  levels = c("Raw", "RPCA")
)

p_lisi <- ggplot(lisi_long_2, aes(x = Method, y = LISI)) +
  geom_boxplot(
    aes(fill = Method),
    width = 0.5,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = c("Raw" = "#E41A1C", "RPCA" = "#377EB8")
  ) +
  labs(
    x = "",
    y = "LISI"
  ) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Raw", "RPCA")
    ),
    label = "p.signif",
    step.increase = -0.12,
    vjust = -0.22
  ) +
  ylim(0, max(lisi_long_2$LISI) * 1.2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

marker_genes <- c(
  # Radial glia
  "PAX6",
  "VIM",
  "GLI3",
  # Endothelial cells
  "CLDN5",
  "PECAM1",
  "VWF",
  "FLT1",
  # Inhibitory neurons
  "GAD1",
  "GAD2",
  "SLC6A1",
  # Oligodendrocyte progenitor cells (OPCs)
  "PDGFRA",
  "CSPG4",
  "OLIG1",
  "OLIG2",
  "SOX10",
  # Microglia
  "CX3CR1",
  "P2RY12",
  "CSF1R",
  # Neuroblasts
  "STMN2",
  # Excitatory neurons
  "SLC17A7",
  "CAMK2A",
  "SATB2",
  # Astrocytes
  "GFAP",
  "AQP4",
  "ALDH1L1",
  "FGFR3",
  "GJA1",
  # Oligodendrocytes
  "MOG",
  "MAG",
  "CLDN11"
)

Idents(objects_plot) <- "CellType"
celltype_levels <- sort(as.character(unique(objects_plot$CellType)))
objects_plot$Celltype <- factor(
  objects_plot$CellType,
  levels = celltype_levels
)

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

p3 <- CellDimPlot(
  objects_plot,
  reduction = "umap.unintegrated",
  group.by = "Dataset",
  palette = "Set1",
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank_axis"
)

p4 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "Dataset",
  palette = "Set1",
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank_axis"
)
p5 <- p3 +
  p4 +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
p6 <- wrap_elements(full = p5) +
  p_lisi +
  plot_layout(widths = c(0.9, 0.1)) +
  plot_annotation(
    tag_levels = "A"
  )
ggsave(
  file.path(fig_dir, "dim_dataset_lisi.pdf"),
  p6,
  width = 12,
  height = 3
)

p7 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "seurat_clusters",
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)

p8 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "Celltype",
  palcolor = color_celltypes[celltype_levels],
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank_axis"
)
p9 <- p7 + p8
ggsave(
  file.path(fig_dir, "dim_seurat_clusters_celltype.pdf"),
  p9,
  width = 15,
  height = 5
)

p13 <- FeatureDimPlot(
  objects_plot,
  features = marker_genes,
  reduction = "umap.rpca",
  ncol = 5,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)
ggsave(
  file.path(fig_dir, "feature_plots.pdf"),
  p13,
  width = 13,
  height = 10
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
  ,
  heatmap_border = TRUE,
  cell_annotation_border = FALSE,
  feature_annotation_border = FALSE,
  heatmap_border_palcolor = "grey",
  heatmap_border_size = 1.5,
  cell_annotation_border_size = 1,
  feature_annotation_border_size = 1,
  height = 8,
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
