#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(scop)
})
source("functions/utils.R")
# Use Cairo for layout measurements as well as final PDF export.
options(device = function(...) grDevices::cairo_pdf(
  file = file.path(tempdir(), "rpca_layout.pdf"), family = "Arial", ...
))

root <- "../../data/BrainOmicsData/integration_25"
validate_annotation_statistics(file.path(root, "heldout_mapping", "rpca_mapping"))
figure_dir <- "figures"

reference <- readRDS(file.path(
  root, "heldout_mapping", "rpca_reference", "rpca_mapping_reference.rds"
))
if (!inherits(reference, "Seurat") || ncol(reference) != 2712452L ||
  !"umap.rpca" %in% Reductions(reference)) {
  stop("Frozen RPCA reference is invalid")
}
reference[["pca"]] <- NULL
reference[["integrated.rpca"]] <- NULL
validate_celltype_metadata(reference[[]], colnames(reference))
reference_levels <- names(brainomics_celltype_colors)[
  names(brainomics_celltype_colors) %in% unique(as.character(reference$CellType))
]
if (!setequal(reference_levels, unique(as.character(reference$CellType)))) {
  stop("Held-out reference CellType labels differ from the manuscript color contract")
}
reference$CellType <- factor(reference$CellType, levels = reference_levels)
gc()

mapping <- fread(
  file.path(
    root, "heldout_mapping", "rpca_mapping", "query_rpca_mapping.tsv.gz"
  ),
  select = c(
    "Cells", "Source_CellType_Display",
    "Ref_UMAP_1", "Ref_UMAP_2"
  )
)
stopifnot(nrow(mapping) == 1655074L, !anyDuplicated(mapping$Cells))
display_mapping <- mapping
display_cells <- nrow(display_mapping)
rm(mapping)
gc()

crosswalk <- fread("integration/celltype_crosswalk.tsv")
source_levels <- crosswalk[
  Label_System == "EGAD00001006049", unique(Display_Label)
]
source_colors <- brainomics_query_celltype_colors
stopifnot(setequal(names(source_colors), source_levels))
source_colors <- source_colors[source_levels]
query_counts <- sparseMatrix(
  i = integer(), j = integer(), dims = c(2L, nrow(display_mapping)),
  dimnames = list(c("placeholderA", "placeholderB"), display_mapping$Cells)
)
query <- CreateSeuratObject(
  query_counts,
  min.cells = 0L,
  min.features = 0L,
  meta.data = data.frame(
    CellType = factor(
      display_mapping$Source_CellType_Display,
      levels = source_levels
    ),
    row.names = display_mapping$Cells
  )
)
query_embedding <- as.matrix(display_mapping[, .(Ref_UMAP_1, Ref_UMAP_2)])
rownames(query_embedding) <- display_mapping$Cells
colnames(query_embedding) <- c("refUMAP_1", "refUMAP_2")
query[["ref.embeddings"]] <- CreateDimReducObject(
  embeddings = query_embedding, key = "refUMAP_", assay = "RNA"
)
rm(display_mapping)
gc()

plot_theme <- list(
  text = ggplot2::element_text(family = "Arial"),
  legend.text = ggplot2::element_text(size = 7.5),
  legend.title = ggplot2::element_text(size = 8.5),
  legend.key.height = grid::unit(3.5, "mm"),
  legend.margin = ggplot2::margin(r = 2, unit = "mm")
)
reference_params <- list(
  palcolor = brainomics_celltype_colors[reference_levels],
  label = FALSE, raster = TRUE, pt.alpha = 1,
  xlab = "UMAP_1", ylab = "UMAP_2",
  theme_use = "theme_blank_axis", theme_args = plot_theme
)
projection_plot <- ProjectionPlot(
  query, reference,
  query_group = "CellType", ref_group = "CellType",
  query_reduction = "ref.embeddings", ref_reduction = "umap.rpca",
  query_param = list(
    palcolor = source_colors, cells.highlight = TRUE, label = FALSE,
    raster = TRUE, xlab = "UMAP_1", ylab = "UMAP_2",
    theme_use = "theme_blank_axis", theme_args = plot_theme
  ),
  ref_param = reference_params,
  pt.size = 0.04, stroke.highlight = 0.03, verbose = FALSE
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  file.path(figure_dir, "rpca_projection.pdf"),
  projection_plot,
  device = grDevices::cairo_pdf, width = 200, height = 85, units = "mm",
  family = "Arial", bg = "white"
)

source("plotting/fig4.R")
