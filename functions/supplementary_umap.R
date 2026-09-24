#!/usr/bin/env Rscript
# Draw one supplementary coverage UMAP from the current 2,602,031-cell atlas.
kind <- commandArgs(trailingOnly = TRUE)
stopifnot(length(kind) == 1L, kind %in% c("age", "region"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
  library(scop)
})
source("functions/config.R")

metadata_path <- file.path(run, "01_metadata", "metadata_working.rds")
embedding_path <- file.path(root, "00_input_audit", "compact", "embedding_umap.rpca.rds")
summary_dir <- file.path(run, "01_metadata")

metadata <- as.data.table(readRDS(metadata_path))
stopifnot(nrow(metadata) == 2602031L, !anyDuplicated(metadata$Cells),
          all(c("Cells", "AgeIntervalID", "BrainRegion") %in% names(metadata)))
metadata <- metadata[, .(Cells, AgeIntervalID, BrainRegion)]
embedding <- readRDS(embedding_path)
stopifnot(is.matrix(embedding), nrow(embedding) == nrow(metadata),
          identical(rownames(embedding), metadata$Cells))

age_levels <- paste0("S", 1:15)
stopifnot(!anyNA(metadata$AgeIntervalID),
          all(metadata$AgeIntervalID %in% age_levels),
          !anyNA(metadata$BrainRegion))
age_counts <- metadata[, .(Cells = .N), by = AgeIntervalID]
region_counts <- metadata[, .(Cells = .N), by = BrainRegion]
expected_age <- fread(file.path(summary_dir, "summary_by_AgeIntervalID.tsv"))
expected_region <- fread(file.path(summary_dir, "summary_by_BrainRegion.tsv"))
stopifnot(nrow(age_counts) == 15L, nrow(region_counts) == 48L,
          all(age_counts$Cells == expected_age$Cells[
            match(age_counts$AgeIntervalID, expected_age$AgeIntervalID)]),
          all(region_counts$Cells == expected_region$Cells[
            match(region_counts$BrainRegion, expected_region$BrainRegion)]))

cells <- metadata$Cells
counts <- sparseMatrix(i = integer(), j = integer(),
                       dims = c(2L, length(cells)),
                       dimnames = list(c("placeholderA", "placeholderB"), cells))
plot_metadata <- data.frame(
  AgeIntervalID = factor(metadata$AgeIntervalID, levels = age_levels),
  BrainRegion = factor(metadata$BrainRegion,
                       levels = sort(unique(metadata$BrainRegion))),
  row.names = cells
)
object <- CreateSeuratObject(
  counts = methods::as(counts, "dgCMatrix"), meta.data = plot_metadata,
  min.cells = 0L, min.features = 0L
)
colnames(embedding) <- c("UMAP_1", "UMAP_2")
object[["umap.rpca"]] <- CreateDimReducObject(
  embeddings = embedding, key = "UMAP_", assay = "RNA"
)
rm(metadata, embedding, counts, plot_metadata)
gc()

age_colors <- setNames(c(
  grDevices::colorRampPalette(c("#0AA344", "#006D87"))(7),
  grDevices::colorRampPalette(c("#2B73AF", "#003D74"))(8)
), age_levels)
region_levels <- levels(object$BrainRegion)
region_index <- seq_along(region_levels)
region_colors <- setNames(grDevices::hcl(
  h = ((region_index - 1L) * 137.508 + 15) %% 360,
  c = rep(c(78, 64, 72), length.out = length(region_levels)),
  l = rep(c(52, 66, 58), length.out = length(region_levels))
), region_levels)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
set.seed(20260906)
if (kind == "age") {
  plot <- scop::CellDimPlot(
    object, reduction = "umap.rpca", group.by = "AgeIntervalID",
    palcolor = age_colors, label = FALSE, raster = TRUE,
    xlab = "UMAP_1", ylab = "UMAP_2", theme_use = "theme_blank_axis"
  )
  output <- "figures/figS1.pdf"
  width <- 125
  height <- 86
} else {
  plot <- scop::CellDimPlot(
    object, reduction = "umap.rpca", group.by = "BrainRegion",
    palcolor = region_colors, label = FALSE, raster = TRUE,
    xlab = "UMAP_1", ylab = "UMAP_2", theme_use = "theme_blank_axis"
  )
  output <- "figures/figS2.pdf"
  plot <- plot & theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.height = grid::unit(2.8, "mm"),
    legend.key.width = grid::unit(3, "mm"),
    legend.margin = margin(2, 0, 0, 0, "mm")
  ) & guides(colour = guide_legend(ncol = 3, byrow = FALSE, title.position = "top"),
             fill = guide_legend(ncol = 3, byrow = FALSE, title.position = "top"))
  width <- 170
  height <- 180
}
ggsave(output, plot, device = grDevices::cairo_pdf,
       width = width, height = height, units = "mm",
       family = "Arial", bg = "white")
message(output, " regenerated from 2,602,031 current cells")
