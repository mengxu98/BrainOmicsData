#!/usr/bin/env Rscript
# Complete literature-defined panel; donor summaries, then equal study weights.
suppressPackageStartupMessages({
  library(scop); library(SeuratObject); library(data.table)
  library(ggplot2); library(patchwork)
})
source("plotting/config.R")
grDevices::pdfFonts(Arial = grDevices::pdfFonts("ArialMT")[[1]])
options(warn = 1)
setDTthreads(2L)
set.seed(20260730)
input <- file.path(doc,"tables/gene_reuse");out<-"figures"
genes <- c("PPP4R2","GXYLT2","KCNJ3","SMCHD1","CTSD","MRPL23","FHIT","SLC10A7","KLHDC4","DACT1","DAAM1")
summary <- fread(file.path(input,"expression_study_equal.tsv"));setnames(summary,"Symbol","Gene")
types <- names(brainomics_celltype_colors)[names(brainomics_celltype_colors)%in%summary$CellType]
stopifnot(setequal(summary$Gene,genes),!anyDuplicated(summary[,.(CellType,Gene)]))
fwrite(summary,file.path(out,"expression_heatmap_source.tsv"),sep="\t")
short <- c("Excitatory neurons" = "Excitatory", "CGE-derived inhibitory neurons" = "CGE inhibitory",
           "MGE-derived inhibitory neurons" = "MGE inhibitory", "Oligodendrocytes" = "Oligodendrocytes",
           "Oligodendrocyte progenitor cells" = "OPC", "Neural progenitor cells" = "Neural progenitor",
           "Oligodendrocyte lineage cells" = "Oligo lineage",
           "Perivascular fibroblasts" = "Perivascular fibroblast",
           "Endothelial cells" = "Endothelial", "Vascular smooth muscle cells" = "Vascular SMC")
display <- setNames(types, types)
display[intersect(names(short), types)] <- short[intersect(names(short), types)]
saveplot <- function(p, name, width = 183, height = 130) {
  for (ext in c("pdf", "png", "svg")) {
    ggsave(file.path(out, paste0(name, ".", ext)), p, width = width, height = height,
           units = "mm", dpi = 300, bg = "white",
           device = switch(ext, pdf = cairo_pdf, png = ragg::agg_png, svg = svglite::svglite))
  }
}
heatmap <- function(value, title, legend, palette) {
  index <- match(paste(rep(types, each = length(genes)), genes), paste(summary$CellType, summary$Gene))
  if (anyNA(index)) stop("Missing gene/type group; must display explicit missingness before plotting")
  mat <- matrix(summary[[value]][index], nrow = length(genes),
                dimnames = list(genes, unname(display[types])))
  if (value == "MeanLogCPM") {
    mat <- t(scale(t(mat)))
    if (any(!is.finite(mat))) stop("A gene has no between-type variation; do not invent a z-score")
  } else mat <- 100 * mat
  color_limits <- if (value == "MeanLogCPM") rep(ceiling(max(abs(mat))), 2) * c(-1, 1) else c(0, 100)
  carrier <- CreateSeuratObject(counts = Matrix::Matrix(0, nrow(mat), ncol(mat),
                                      sparse = TRUE, dimnames = dimnames(mat)),
                                assay = "Summary",
                                meta.data = data.frame(Type = factor(colnames(mat), levels = colnames(mat)),
                                  Unit = "study_equal_donor_summary", row.names = colnames(mat)),
                                min.cells = 0, min.features = 0)
  LayerData(carrier, assay = "Summary", layer = "data") <- mat
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  h <- scop::FeatureHeatmap(carrier, features = genes, group.by = NULL, assay = "Summary",
    max_cells = ncol(mat), cell_order = colnames(mat),
    layer = "data", lib_normalize = FALSE, exp_method = "raw", exp_legend_title = legend,
    heatmap_palette = palette, limits = color_limits, border = FALSE,
    nlabel = 0, show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left",
    row_names_rot = 0, column_names_rot = 45, column_names_side = "bottom", column_title = title,
    cluster_rows = FALSE, cluster_columns = FALSE, anno_terms = FALSE,
    anno_keys = FALSE, anno_features = FALSE, use_raster = FALSE,
    width = 2.25, height = 1.55, legend.position = "bottom", verbose = FALSE,
    ht_params = list(row_labels = ifelse(genes == "CTSD", "CTSD \u2020", genes),
      row_names_gp = grid::gpar(fontfamily = "Arial", fontsize = 7, fontface = "italic"),
      column_names_gp = grid::gpar(fontfamily = "Arial", fontsize = 6.5),
      column_names_max_height = grid::unit(40, "mm"),
      column_title_gp = grid::gpar(fontfamily = "Arial", fontsize = 8)))
  observed <- h$matrix_list[[1]][rownames(mat), colnames(mat), drop = FALSE]
  stopifnot(isTRUE(all.equal(as.numeric(observed), as.numeric(mat), tolerance = 1e-7)))
  h$plot
}
a <- heatmap("MeanLogCPM", "Expression across cell types", "Expression z-score", "RdBu")
b <- heatmap("Detection", "Within-type detection", "Detected cells (%)", "viridis")
caption <- "S15 (age 60+), prefrontal cortex. Donor means within study, then equal study weights.\nA: row z-score of mean log1p CPM. All 11 Table S6 candidates; absolute expression and coverage in source data.\n\u2020 CTSD: very low detection in ROSMAP; relative color intensity does not imply high expression."
saveplot((a | b) + plot_annotation(caption = caption, tag_levels = "A",
  theme = theme(text = element_text(family = "Arial", size = 8),
                plot.caption = element_text(size = 7, hjust = 0))), "reuse_gene_expression", 210, 130)

study<-fread(file.path(input,"effects_by_study.tsv"));pooled<-fread(file.path(input,"effects_study_equal.tsv"))
setnames(study,c("Symbol","Effect"),c("Gene","Estimate"));setnames(pooled,c("Symbol","Effect"),c("Gene","Estimate"));pooled[,Dataset:="Study-equal"]
effects <- rbindlist(list(study, pooled), fill = TRUE)
stopifnot(setequal(unique(effects$Gene), genes), !anyDuplicated(effects[, .(Dataset, Gene)]))
effects[, Gene := factor(Gene, levels = rev(genes))]
dataset_levels <- c(sort(unique(study$Dataset)), "Study-equal")
effects[, Dataset := factor(Dataset, levels = dataset_levels)]
fwrite(effects, file.path(out, "oligo_microglia_effect_source.tsv"), sep = "\t")
finite <- is.finite(effects$Estimate)
if (!all(finite)) message(sum(!finite), " effects not estimable; listed in source table")
pal <- setNames(c("#4575B4", "#D68A43", "#262626"), dataset_levels)
dodge <- position_dodge(width = .65)
p <- ggplot(effects[finite], aes(x = Estimate, y = Gene, colour = Dataset, group = Dataset)) +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey65", linetype = "dashed") +
  geom_errorbar(aes(xmin = CI95_Lower, xmax = CI95_Upper), orientation = "y",
                width = .14, linewidth = .4, position = dodge, na.rm = TRUE) +
  geom_point(aes(shape = Dataset), size = 1.5, position = dodge) +
  scale_colour_manual(values = pal) + scale_shape_manual(values = c(16, 17, 18)) +
  scale_y_discrete(drop = FALSE, labels = function(z) ifelse(z == "CTSD", "CTSD \u2020", z)) +
  labs(x = "Oligodendrocytes minus microglia (log1p CPM)", y = NULL,
       title = "Cell-type contrasts across studies", colour = NULL, shape = NULL,
       caption = "All 11 candidates; 40 paired donors in two studies. Bars: pointwise t intervals across the two study means.\nNo resampling; study-specific estimates have no uncertainty bars. Full coverage and contrasts are in source data.\n\u2020 CTSD: very low detection in ROSMAP; zero contrasts do not establish equivalence.") +
  scop::theme_scop() + theme(text = element_text(family = "Arial", size = 8),
     axis.text.y = element_text(face = "italic", size = 8),
     axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 8),
     plot.title = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 7),
     plot.caption = element_text(size = 7, hjust = 0))
saveplot(p, "reuse_gene_study_effects", 183, 160)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
