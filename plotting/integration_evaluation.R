suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(Seurat)
  library(SeuratObject)
  library(scop)
})

source("functions/utils.R")

evaluation_dir <- "../../data/BrainOmicsData/integration_25/evaluation"
lisi_input <- file.path(evaluation_dir, "dataset_lisi_dataset_equal.tsv")
figure_dir <- "figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

dataset_scores <- read_tsv(lisi_input)

expected_columns <- c("Dataset", "Known_Donors", method_levels)
if (!all(expected_columns %in% names(dataset_scores))) {
  stop("The LISI source table does not contain the required columns")
}
if (nrow(dataset_scores) != 25L || length(unique(dataset_scores$Dataset)) != 25L) {
  stop("The revised LISI panel requires exactly 25 independent datasets")
}
if (anyNA(dataset_scores[, method_levels, drop = FALSE]) ||
  any(!is.finite(as.matrix(dataset_scores[, method_levels, drop = FALSE])))) {
  stop("The four-method LISI source contains missing or non-finite values")
}

lisi_long <- do.call(rbind, lapply(method_levels, function(method) {
  data.frame(
    Dataset = dataset_scores$Dataset,
    Known_Donors = dataset_scores$Known_Donors,
    Method = method,
    LISI = dataset_scores[[method]],
    stringsAsFactors = FALSE
  )
}))
lisi_long$Method <- factor(lisi_long$Method, levels = method_levels)

theme_set(
  theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5.5),
      strip.text = element_text(size = 6.5, face = "bold"),
      plot.title = element_text(size = 7, face = "bold"),
      plot.tag = element_text(size = 8, face = "bold")
    )
)

p_lisi <- ggplot(lisi_long, aes(x = Method, y = LISI)) +
  geom_boxplot(
    aes(fill = Method),
    width = 0.56, outlier.shape = NA,
    linewidth = 0.36, colour = "#303030"
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.07))
  ) +
  labs(
    x = NULL,
    y = "iLISI"
  ) +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.36, colour = "#303030"),
    axis.text = element_text(colour = "#303030"),
    axis.text.x = element_text(size = 6.2),
    axis.title.y = element_text(size = 6.6),
    plot.margin = margin(2, 2, 1, 2)
  )
ggplot2::ggsave(
  file.path(figure_dir, "fig2_lisi.pdf"), p_lisi,
  device = grDevices::cairo_pdf, width = 62, height = 66, units = "mm",
  family = "Arial", bg = "white"
)
write_tsv(lisi_long, file.path(figure_dir, "fig2_lisi_data.tsv"))

rm(p_lisi)
gc()

umap <- readRDS(file.path(evaluation_dir, "umap_plot_data.rds"))
dataset_levels <- levels(umap$Dataset)
if (is.null(dataset_levels)) {
  dataset_levels <- sort(unique(as.character(umap$Dataset)))
  umap$Dataset <- factor(umap$Dataset, levels = dataset_levels)
}
# Keep the complete horizontal panel legible when placed in the main figure.
umap_panel_theme <- theme(
  plot.title = element_text(size = 10, face = "bold"),
  plot.subtitle = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10)
)

if (nrow(umap) != 2712452L || anyDuplicated(umap$Cell) ||
  anyNA(umap[, paste0(rep(method_levels, each = 2L), "_", 1:2)])) {
  stop("The complete four-method UMAP table violates the frozen cell contract")
}
placeholder <- methods::as(Matrix::sparseMatrix(
  i = integer(), j = integer(), x = numeric(),
  dims = c(2L, nrow(umap)),
  dimnames = list(c("placeholderA", "placeholderB"), umap$Cell)
), "dgCMatrix")
umap_object <- CreateSeuratObject(
  counts = placeholder,
  meta.data = data.frame(Dataset = umap$Dataset, row.names = umap$Cell),
  min.cells = 0L,
  min.features = 0L
)
embedding <- as.matrix(umap[, c("Raw_1", "Raw_2")])
rownames(embedding) <- umap$Cell
colnames(embedding) <- paste0("UMAP_", 1:2)
umap_object[["plot.umap"]] <- CreateDimReducObject(
  embeddings = embedding, key = "UMAP_", assay = "RNA"
)
p_umap_raw <- scop::CellDimPlot(
  umap_object,
  reduction = "plot.umap",
  group.by = "Dataset",
  palcolor = dataset_colors[dataset_levels],
  label = FALSE,
  raster = TRUE,
  raster.dpi = c(1600, 1600),
  show_stat = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  title = "Raw",
  legend.position = "none",
  theme_use = theme_blank_axis,
  theme_args = list(lab_size = 10, axis_lwd = 1)
) +
  umap_panel_theme
rm(embedding)
gc()

embedding <- as.matrix(umap[, c("scVI_1", "scVI_2")])
rownames(embedding) <- umap$Cell
colnames(embedding) <- paste0("UMAP_", 1:2)
umap_object[["plot.umap"]] <- CreateDimReducObject(
  embeddings = embedding, key = "UMAP_", assay = "RNA"
)
p_umap_scvi <- scop::CellDimPlot(
  umap_object,
  reduction = "plot.umap",
  group.by = "Dataset",
  palcolor = dataset_colors[dataset_levels],
  label = FALSE,
  raster = TRUE,
  raster.dpi = c(1600, 1600),
  show_stat = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  title = "scVI",
  legend.position = "none",
  theme_use = theme_blank_axis,
  theme_args = list(lab_size = 10, axis_lwd = 1)
) +
  umap_panel_theme

rm(embedding)
gc()

embedding <- as.matrix(umap[, c("Harmony_1", "Harmony_2")])
rownames(embedding) <- umap$Cell
colnames(embedding) <- paste0("UMAP_", 1:2)
umap_object[["plot.umap"]] <- CreateDimReducObject(
  embeddings = embedding, key = "UMAP_", assay = "RNA"
)
p_umap_harmony <- scop::CellDimPlot(
  umap_object,
  reduction = "plot.umap",
  group.by = "Dataset",
  palcolor = dataset_colors[dataset_levels],
  label = FALSE,
  raster = TRUE,
  raster.dpi = c(1600, 1600),
  show_stat = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  title = "Harmony",
  legend.position = "none",
  theme_use = theme_blank_axis,
  theme_args = list(lab_size = 10, axis_lwd = 1)
) +
  umap_panel_theme

rm(embedding)
gc()

embedding <- as.matrix(umap[, c("RPCA_1", "RPCA_2")])
rownames(embedding) <- umap$Cell
colnames(embedding) <- paste0("UMAP_", 1:2)
umap_object[["plot.umap"]] <- CreateDimReducObject(
  embeddings = embedding, key = "UMAP_", assay = "RNA"
)
p_umap_rpca <- scop::CellDimPlot(
  umap_object,
  reduction = "plot.umap",
  group.by = "Dataset",
  palcolor = dataset_colors[dataset_levels],
  label = FALSE,
  raster = TRUE,
  raster.dpi = c(1600, 1600),
  show_stat = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  title = "RPCA",
  legend.position = "right",
  theme_use = theme_blank_axis,
  theme_args = list(lab_size = 10, axis_lwd = 1)
) +
  umap_panel_theme

rm(embedding)
gc()
rm(placeholder, umap_object, umap)
gc()

p_umap_comparison <- patchwork::wrap_plots(
  p_umap_raw, p_umap_scvi, p_umap_harmony, p_umap_rpca,
  nrow = 1, widths = c(1, 1, 1, 1), guides = "keep"
)
ggplot2::ggsave(
  file.path(figure_dir, "integration_umap.pdf"), p_umap_comparison,
  device = grDevices::cairo_pdf, width = 312, height = 70, units = "mm",
  family = "Arial", bg = "white"
)

rm(p_umap_raw, p_umap_scvi, p_umap_harmony, p_umap_rpca, p_umap_comparison)
gc()
message("Four-method LISI panel and integration evaluation written to ", figure_dir)
