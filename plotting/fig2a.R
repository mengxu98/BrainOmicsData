#!/usr/bin/env Rscript
# Figure 2A: draw saved UMAP coordinates without fitting any embedding.
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(scop)
})
source("functions/config.R")

metadata <- fread(resolve_input_file(
  "BRAINOMICS_FIG2_METADATA_FILE", root, "core_metadata_minimal.tsv.gz"
))
stopifnot(nrow(metadata) == 2602031L, !anyDuplicated(metadata$Cells),
          uniqueN(metadata$Dataset) == 22L, !anyNA(metadata$Dataset))
datasets <- sort(unique(metadata$Dataset))
stopifnot(all(datasets %in% names(dataset_colors)))
reductions <- c(Raw = "umap.unintegrated", scVI = "umap.scvi",
                Harmony = "umap.harmony", RPCA = "umap.rpca")
variables <- c(Raw = "BRAINOMICS_RAW_UMAP_FILE", scVI = "BRAINOMICS_SCVI_UMAP_FILE",
               Harmony = "BRAINOMICS_HARMONY_UMAP_FILE", RPCA = "BRAINOMICS_RPCA_UMAP_FILE")
paths <- vapply(method_levels, function(method) resolve_input_file(
  variables[[method]], root, paste0("embedding_", reductions[[method]], ".rds")
), character(1))

cells <- metadata$Cells
counts <- sparseMatrix(i = integer(), j = integer(), x = numeric(),
  dims = c(2L, length(cells)),
  dimnames = list(c("placeholderA", "placeholderB"), cells))
object <- CreateSeuratObject(counts = counts,
  meta.data = data.frame(Dataset = factor(metadata$Dataset, levels = datasets),
                        row.names = cells), min.cells = 0L, min.features = 0L)
set.seed(20260906)
draw_order <- cells[sample.int(length(cells))]
plots <- vector("list", length(method_levels))
for (i in seq_along(method_levels)) {
  method <- method_levels[[i]]
  embedding <- readRDS(paths[[method]])
  stopifnot(is.matrix(embedding), ncol(embedding) == 2L,
            identical(rownames(embedding), cells), all(is.finite(embedding)))
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  object[["display"]] <- CreateDimReducObject(
    embeddings = embedding, key = "UMAP_", assay = "RNA")
  p <- scop::CellDimPlot(object, reduction = "display", group.by = "Dataset",
    palcolor = dataset_colors[datasets], label = FALSE, seed = 11,
    raster = TRUE, raster.dpi = c(1600, 1600), pt.size = 2,
    xlab = "UMAP_1", ylab = "UMAP_2", legend.title = "Dataset",
    theme_use = "theme_blank_axis", combine = FALSE)[[1]]
  p <- order_dim_plot_cells(p, draw_order)
  stopifnot(identical(rownames(p$data), draw_order))
  p$layers <- Filter(function(layer) !inherits(layer$geom, "GeomCustomAnn"), p$layers)
  plots[[i]] <- p + theme_blank_axis(lab_size = 5, axis_lwd = .6) +
    p$coordinates + labs(title = method, subtitle = paste0("nCells:", length(cells))) +
    theme(text = element_text(family = "Arial", size = 6),
      plot.title = element_text(size = 6, face = "bold"),
      plot.subtitle = element_text(size = 5),
      legend.text = element_text(size = 5), legend.title = element_text(size = 5.5),
      legend.key.height = grid::unit(2.3, "mm"),
      legend.key.width = grid::unit(2, "mm"),
      plot.margin = margin(1, 1, 1, 1, unit = "mm")) +
    guides(colour = guide_legend(ncol = 2, byrow = FALSE,
      override.aes = list(size = 1.1)))
  message("Figure 2A coordinates checked: ", method, " / ", length(cells), " cells")
}
panel <- wrap_plots(c(plots, list(guide_area())), nrow = 1,
                   widths = c(1, 1, 1, 1, 2.1), guides = "collect")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/fig2a.pdf", panel, device = cairo_pdf,
       width = 176, height = 37, units = "mm", family = "Arial", bg = "white")
