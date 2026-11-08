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
  plots[[i]] <- p + theme_void() + labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  message("Figure 2A coordinates checked: ", method, " / ", length(cells), " cells")
}
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
grDevices::cairo_pdf("figures/fig2a.pdf", width = 176 / 25.4, height = 68 / 25.4, family = "Arial")
grid::grid.newpage()
for (i in seq_along(method_levels)) {
  x <- (i - .5) / 4
  print(plots[[i]], newpage = FALSE, vp = grid::viewport(x = x, y = grid::unit(42.5, "mm"),
    width = grid::unit(39, "mm"), height = grid::unit(39 * 63.424983 / 74.328693, "mm")))
  grid::grid.text(method_levels[i], x = x, y = grid::unit(65, "mm"),
    gp = grid::gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold"))
  x0 <- (i - 1) * 44 + 2
  grid::grid.lines(grid::unit(c(x0, x0 + 5), "mm"), grid::unit(c(26.5, 26.5), "mm"),
    arrow = grid::arrow(length = grid::unit(1, "mm")), gp = grid::gpar(lwd = .6))
  grid::grid.lines(grid::unit(c(x0, x0), "mm"), grid::unit(c(26.5, 31.5), "mm"),
    arrow = grid::arrow(length = grid::unit(1, "mm")), gp = grid::gpar(lwd = .6))
  grid::grid.text("UMAP_1", grid::unit(x0, "mm"), grid::unit(24.5, "mm"), just = "left",
    gp = grid::gpar(fontfamily = "Arial", fontsize = 5.8))
  grid::grid.text("UMAP_2", grid::unit(x0 - 1.5, "mm"), grid::unit(28.5, "mm"), rot = 90,
    gp = grid::gpar(fontfamily = "Arial", fontsize = 5.8))
}
grid::grid.text(paste0("n = ", format(length(cells), big.mark = ",", scientific = FALSE), " cells"),
  x = .125, y = grid::unit(60.5, "mm"), gp = grid::gpar(fontfamily = "Arial", fontsize = 6.4, col = "#444444"))
grid::grid.text("Dataset (cells)", grid::unit(2, "mm"), grid::unit(18.5, "mm"), just = "left",
  gp = grid::gpar(fontfamily = "Arial", fontsize = 6.2, fontface = "bold"))
counts <- table(factor(metadata$Dataset, levels = datasets))
for (i in seq_along(datasets)) {
  x <- (i %/% 6) * 44 + 2
  y <- 18.5 - (i %% 6) * 3.2
  grid::grid.points(grid::unit(x, "mm"), grid::unit(y, "mm"), pch = 16,
    size = grid::unit(1.6, "mm"), gp = grid::gpar(col = dataset_colors[datasets[i]]))
  label <- paste0(datasets[i], " (", format(counts[i], big.mark = ",", scientific = FALSE, trim = TRUE), ")")
  grid::grid.text(label, grid::unit(x + 2, "mm"), grid::unit(y, "mm"), just = "left",
    gp = grid::gpar(fontfamily = "Arial", fontsize = 5.9))
}
grDevices::dev.off()
