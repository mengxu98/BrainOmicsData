#!/usr/bin/env Rscript
# Original-source labels on existing UMAPs; raw refers to label origin.
# External and independent supplementary annotations remain separate from this plot.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
  library(scop)
})
source("functions/utils.R")
root <- "../../data/BrainOmicsData/integration_25"
umap <- readRDS(file.path(root, "evaluation/umap_plot_data.rds"))
labels <- as.data.table(readRDS(file.path(root, "annotation/source_labels_harmonized.rds")))
cells <- as.character(umap$Cell)
index <- match(cells, labels$Cells)
stopifnot(length(cells) == 2712452L, !anyDuplicated(cells),
          !anyDuplicated(labels$Cells), nrow(labels) == length(cells), !anyNA(index))
types <- labels$Harmonized_CellType[index]
unavailable <- labels$Mapping_Status[index] == "source_label_unavailable"
stopifnot(!anyNA(unavailable), sum(unavailable) == 53000L,
          all(is.na(types[unavailable])))
# Display-only grouping; original annotation and mapping status remain unchanged.
unassigned <- unavailable | types %in% c("Source-reported unknown", "Mixed or low-quality")
stopifnot(sum(unassigned) == 74745L)
types[unassigned] <- "Source annotation unavailable"
types[types == "Source-reported Immune Myocytes"] <- "Immune Myocytes"
stopifnot(!anyNA(types), all(nzchar(types)))
type_order <- sort(unique(types))
neutral <- c("Immune Myocytes", "Source annotation unavailable")
type_order <- c(setdiff(type_order, neutral), intersect(neutral, type_order))
colors <- setNames(grDevices::hcl.colors(length(type_order), "Dark 3"), type_order)
# Align identical biological identities with the existing cell-type panel.
shared <- c("Astrocyte" = "Astrocytes", "Endothelial cell" = "Endothelial cells",
            "Excitatory neuron" = "Excitatory neurons", "Microglia" = "Microglia",
            "Oligodendrocyte" = "Oligodendrocytes", "Radial glia" = "Radial glia",
            "Oligodendrocyte precursor cell" = "Oligodendrocyte progenitor cells")
colors[names(shared)] <- brainomics_celltype_colors[shared]
colors["Inhibitory neuron"] <- "#286A9D"
colors["Immune Myocytes"] <- "#4E4E4E"
colors["Source annotation unavailable"] <- "#C7C7C7"

counts <- sparseMatrix(i = integer(), j = integer(), x = numeric(),
                       dims = c(2L, length(cells)),
                       dimnames = list(c("placeholderA", "placeholderB"), cells))
object <- CreateSeuratObject(counts = counts,
  meta.data = data.frame(SourceType = factor(types, levels = type_order), row.names = cells),
  min.cells = 0L, min.features = 0L)
set.seed(42)
draw_order <- cells[sample.int(length(cells))]
methods <- c("Raw", "scVI", "Harmony", "RPCA")
plots <- lapply(methods, function(method) {
  embedding <- as.matrix(umap[, paste0(method, c("_1", "_2"))])
  stopifnot(all(is.finite(embedding)))
  dimnames(embedding) <- list(cells, c("UMAP_1", "UMAP_2"))
  object[["umap.compare"]] <- CreateDimReducObject(
    embeddings = embedding, key = "COMPARE_", assay = "RNA")
  p <- scop::CellDimPlot(object, reduction = "umap.compare", group.by = "SourceType",
    palcolor = colors, label = FALSE, seed = 11,
    raster = TRUE, raster.dpi = c(1600, 1600), pt.size = 2,
    title = method, xlab = "UMAP_1", ylab = "UMAP_2",
    legend.title = "Source cell type", legend.position = "bottom",
    theme_use = "theme_blank_axis",
    theme_args = list(text = element_text(family = "Arial")), combine = FALSE)[[1]]
  p <- order_dim_plot_cells(p, draw_order)
  stopifnot(identical(rownames(p$data), draw_order), nrow(p$data) == length(cells))
  p + theme(plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 10),
            legend.text = element_text(size = 7), legend.title = element_text(size = 9)) +
    guides(colour = guide_legend(ncol = 4, byrow = TRUE))
})
names(plots) <- methods
plot <- patchwork::wrap_plots(plots[c("scVI", "Harmony", "RPCA")], nrow = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.box.margin = margin(t = 8, unit = "mm"))
ggsave("figures/integration_celltype_umap_raw.pdf", plot,
       device = grDevices::cairo_pdf, width = 312, height = 154, units = "mm",
       family = "Arial", bg = "white")
plot_four <- patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.box.margin = margin(t = 8, unit = "mm"))
ggsave("figures/integration_celltype_umap_raw_4methods.pdf", plot_four,
       device = grDevices::cairo_pdf, width = 390, height = 154, units = "mm",
       family = "Arial", bg = "white")
message("All ", length(cells), " cells displayed; ", sum(unassigned),
        " cells shown as Source annotation unavailable (53,000 unavailable, 19,187 unknown, 2,558 mixed/low-quality); no supplementary labels used.")
