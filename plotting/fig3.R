#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scop)
  library(Seurat)
})

source("functions/data_paths.R")
source("functions/metadata_schema.R")
source("functions/dataset_metadata.R")
source("functions/utils.R")
validate_annotation_statistics(
  "../../data/BrainOmicsData/integration_25/evaluation/reference_summary"
)

result_root <- "../../data/BrainOmicsData/integration_25"
figure_dir <- "figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
summary_dir <- file.path(result_root, "evaluation", "reference_summary")
age <- fread(file.path(summary_dir, "age_interval_summary.tsv"))
coverage <- fread(file.path(summary_dir, "celltype_age_region_coverage.tsv"))
assignments <- as.data.table(read_celltype_assignments())
plot_metadata <- apply_curated_age_interval_decisions(
  readRDS(file.path(result_root, "evaluation", "plot_metadata_slim.rds"))
)
if (
  nrow(assignments) != 2712452L ||
    anyDuplicated(assignments$Cells) ||
    length(unique(assignments$CellType)) != 16L
) {
  stop("Figure 3 assignments differ from the frozen 16-class contract")
}

message("[fig3] Using all-cell UMAP coordinates without loading expression")
umap <- readRDS(file.path(result_root, "evaluation", "umap_plot_data.rds"))
cells <- as.character(umap$Cell)
stopifnot(
  identical(cells, assignments$Cells),
  identical(cells, as.character(plot_metadata$Cells)),
  !anyDuplicated(cells)
)
object <- CreateSeuratObject(
  counts = Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    dims = c(2L, length(cells)),
    dimnames = list(c("placeholderA", "placeholderB"), cells)
  ),
  min.cells = 0L,
  min.features = 0L
)
coordinates <- as.matrix(umap[, c("RPCA_1", "RPCA_2")])
dimnames(coordinates) <- list(cells, c("UMAP_1", "UMAP_2"))
object[["umap.rpca"]] <- CreateDimReducObject(
  embeddings = coordinates,
  key = "UMAP_",
  assay = "RNA"
)
rm(umap, coordinates)
gc()
object$CellType <- assignments$CellType
object$AgeIntervalID <- as.character(plot_metadata$AgeIntervalID)
object$BrainRegion <- standardize_brain_region(plot_metadata$BrainRegion)
age_levels <- paste0("S", 1:15)
unknown_age <- is.na(object$AgeIntervalID) |
  !object$AgeIntervalID %in% age_levels
if (any(unknown_age)) {
  stop("Figure 3 metadata contains unresolved age intervals")
}
object$AgeIntervalID <- factor(object$AgeIntervalID, levels = age_levels)

age_colors <- setNames(
  c(
    grDevices::colorRampPalette(c("#0AA344", "#006D87"))(7),
    grDevices::colorRampPalette(c("#2B73AF", "#003D74"))(8)
  ),
  age_levels
)
p_age_umap <- scop::CellDimPlot(
  object,
  reduction = "umap.rpca",
  group.by = "AgeIntervalID",
  palcolor = age_colors,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank_axis"
)
ggplot2::ggsave(
  file.path(figure_dir, "fig3_age_umap.pdf"),
  p_age_umap,
  device = grDevices::cairo_pdf,
  width = 125,
  height = 86,
  units = "mm",
  family = "Arial",
  bg = "white"
)

region_levels <- sort(unique(as.character(object$BrainRegion)))
region_index <- seq_along(region_levels)
region_colors <- setNames(
  grDevices::hcl(
    h = ((region_index - 1L) * 137.508 + 15) %% 360,
    c = rep(c(78, 64, 72), length.out = length(region_levels)),
    l = rep(c(52, 66, 58), length.out = length(region_levels))
  ),
  region_levels
)
object$BrainRegion <- factor(object$BrainRegion, levels = region_levels)
p_region_umap <- scop::CellDimPlot(
  object,
  reduction = "umap.rpca",
  group.by = "BrainRegion",
  palcolor = region_colors,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank_axis"
)
ggplot2::ggsave(
  file.path(figure_dir, "fig3_region_umap.pdf"),
  p_region_umap,
  device = grDevices::cairo_pdf,
  width = 330,
  height = 100,
  units = "mm",
  family = "Arial",
  bg = "white"
)

age <- age[AgeIntervalID %in% age_levels]

age$AgeIntervalID <- factor(age$AgeIntervalID, levels = age_levels)
setorder(age, AgeIntervalID)
age_table <- age[, .(
  `Age interval ID` = as.character(AgeIntervalID),
  `Age interval` = AgeInterval,
  `Boundary` = AgeRange,
  `Cells/nuclei` = format(Cells, big.mark = ",", scientific = FALSE),
  `Known donors` = Known_Donors,
  `Datasets` = Datasets,
  `Brain regions` = Brain_Regions
)]
age_table <- as.data.frame(age_table, stringsAsFactors = FALSE)
table_text_gp <- grid::gpar(fontsize = 8)
table_padding_mm <- 1.5
column_widths_mm <- vapply(
  seq_len(ncol(age_table)),
  function(j) {
    labels <- c(names(age_table)[j], as.character(age_table[[j]]))
    widths <- lapply(
      labels,
      function(label) grid::grobWidth(grid::textGrob(label, gp = table_text_gp))
    )
    grid::convertWidth(max(do.call(grid::unit.c, widths)), "mm", TRUE) +
      table_padding_mm * 2
  },
  numeric(1)
)
color_width_mm <- 5
row_height_mm <- 5.5
header_height_mm <- 5
table_width_mm <- sum(column_widths_mm)
plot_width_mm <- color_width_mm + table_width_mm
plot_height_mm <- row_height_mm * nrow(age_table) + header_height_mm

draw_age_table <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = grid::unit(2, "mm"),
    y = grid::unit(1, "mm"),
    width = grid::unit(plot_width_mm, "mm"),
    height = grid::unit(plot_height_mm, "mm"),
    just = c("left", "bottom")
  ))

  n_rows <- nrow(age_table)
  grid::grid.rect(
    x = grid::unit(0, "mm"),
    y = grid::unit(
      header_height_mm + row_height_mm * (n_rows - seq_len(n_rows)),
      "mm"
    ),
    width = grid::unit(color_width_mm, "mm"),
    height = grid::unit(row_height_mm, "mm"),
    gp = grid::gpar(
      fill = age_colors[as.character(age_table[["Age interval ID"]])],
      col = NA
    ),
    just = c("left", "bottom")
  )

  column_left_mm <- color_width_mm
  for (j in seq_len(ncol(age_table))) {
    for (i in seq_len(n_rows)) {
      grid::grid.rect(
        x = grid::unit(column_left_mm, "mm"),
        y = grid::unit(
          header_height_mm + row_height_mm * (n_rows - i),
          "mm"
        ),
        width = grid::unit(column_widths_mm[j], "mm"),
        height = grid::unit(row_height_mm, "mm"),
        gp = grid::gpar(fill = "#EFF7FC", col = "#9FBBd3", lwd = 0.5),
        just = c("left", "bottom")
      )
      grid::grid.text(
        as.character(age_table[i, j]),
        x = grid::unit(column_left_mm + table_padding_mm, "mm"),
        y = grid::unit(
          header_height_mm + row_height_mm * (n_rows - i + 0.5),
          "mm"
        ),
        gp = table_text_gp,
        just = "left"
      )
    }
    grid::grid.text(
      names(age_table)[j],
      x = grid::unit(column_left_mm + column_widths_mm[j] / 2, "mm"),
      y = grid::unit(header_height_mm / 2, "mm"),
      gp = table_text_gp
    )
    column_left_mm <- column_left_mm + column_widths_mm[j]
  }
  grid::popViewport()
}
grDevices::cairo_pdf(
  file.path(figure_dir, "fig3_age_table.pdf"),
  width = (plot_width_mm + 4) / 25.4,
  height = (plot_height_mm + 2) / 25.4,
  family = "Arial",
  onefile = FALSE
)
draw_age_table()
grDevices::dev.off()
age_celltype <- coverage[
  AgeIntervalID %in% age_levels,
  .(Cells = sum(Cells)),
  by = .(CellType, AgeIntervalID)
]
age_celltype <- merge(
  CJ(CellType = unique(coverage$CellType), AgeIntervalID = age_levels),
  age_celltype,
  by = c("CellType", "AgeIntervalID"),
  all.x = TRUE
)
age_celltype[is.na(Cells), Cells := 0L]
# The coverage table is already donor-aware.  The heatmap uses cell abundance
# only for descriptive coverage and does not perform inferential statistics.
stopifnot(all(is.finite(age_celltype$Cells)), all(age_celltype$Cells >= 0))
age_celltype[, Log10_Cells := log10(pmax(Cells, 0) + 1)]
age_celltype[, AgeIntervalID := factor(AgeIntervalID, levels = age_levels)]
p_coverage <- ggplot(
  age_celltype,
  aes(x = AgeIntervalID, y = CellType, fill = Log10_Cells)
) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours = c("#EFF7FC", "#0AA344", "#006D87", "#003D74"),
    name = expression(log[10](cells + 1))
  ) +
  labs(
    x = "Age interval ID",
    y = "Cell type",
    title = "Cell-type coverage across the 15 age intervals"
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6.0),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold")
  )
ggplot2::ggsave(
  file.path(figure_dir, "fig3_celltype_age.pdf"),
  p_coverage,
  device = grDevices::cairo_pdf,
  width = 115,
  height = 88,
  units = "mm",
  family = "Arial",
  bg = "white"
)

message("Figure 3 individual panels completed in ", figure_dir)
