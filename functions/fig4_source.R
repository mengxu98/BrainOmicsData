#!/usr/bin/env Rscript
# Figure 4 source mapping and donor concordance panels.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrastr);library(ComplexHeatmap)})
source("functions/config.R")
grDevices::pdfFonts(Arial=grDevices::pdfFonts('ArialMT')[[1]])
options(device=function(...) grDevices::cairo_pdf(file=file.path(tempdir(),'egad-layout.pdf'),family='Arial',...))
results_root<-file.path(doc,'tables/egad_mapping')
required_inputs <- file.path(results_root, c(
  'query_predictions.tsv.gz',
  'full_source_prediction_counts.tsv',
  'comparable_donor_source_agreement.tsv'
))
stopifnot(all(file.exists(required_inputs)))
out<-'figures';dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-fread(file.path(results_root,'query_predictions.tsv.gz'))
setnames(d,c('CellID','Source_CellClass','Working_CellType'),c('Cells','Source_Label','Predicted_CellType'))
stopifnot(nrow(d)==1661798L,!anyDuplicated(d$Cells),all(d$Mapping_Status=='mapped'),all(is.finite(d$UMAP_1)),all(is.finite(d$UMAP_2)))
stopifnot(all(d$Predicted_CellType %in% names(brainomics_celltype_colors)),uniqueN(d$Predicted_CellType)==12L)
# Source names are display-only harmonizations; coordinates and votes are unchanged.
rename<-c('Neuron'='Neurons','Neuroblast'='Neuroblasts','Neuronal IPC'='Neuronal intermediate progenitor cells','Glioblast'='Glioblasts','Oligo'='Oligodendrocyte lineage','Immune'='Immune cells','Vascular'='Vascular cells','Fibroblast'='Perivascular fibroblasts','Erythrocyte'='Erythrocytes')
d[Source_Label%in%names(rename),Source_Label:=unname(rename[Source_Label])]
# Reuse Fig. 2G's CellDimPlot rendering and final-size arrow-axis theme.
suppressPackageStartupMessages({library(SeuratObject);library(scop)})
qcolors<-brainomics_query_celltype_colors
rcolors<-c(brainomics_celltype_colors,'Not mapped (no counts)'='#999999',
 'Not mapped (no reference-feature counts)'='#999999')
counts<-Matrix::sparseMatrix(i=integer(),j=integer(),dims=c(2L,nrow(d)),
 dimnames=list(c('displayA','displayB'),d$Cells))
metadata<-data.frame(Source_Label=factor(d$Source_Label,levels=names(qcolors)),
 Predicted_CellType=factor(d$Predicted_CellType,levels=names(rcolors)),row.names=d$Cells)
metadata<-droplevels(metadata)
object<-CreateSeuratObject(counts=counts,meta.data=metadata)
embedding<-as.matrix(d[,.(UMAP_1,UMAP_2)]);rownames(embedding)<-d$Cells
object[['display']]<-CreateDimReducObject(embeddings=embedding,key='UMAP_',assay='RNA')
set.seed(2026);draw_order<-d$Cells[sample.int(nrow(d))]
# Keep 1:1 axis units, and size the cloud to the cell-type legend.
umap_xlim <- range(d$UMAP_1); umap_ylim <- range(d$UMAP_2)
aspect_ratio <- diff(umap_xlim) / diff(umap_ylim)
umap_panel_height <- 34
umap_panel_width <- round(umap_panel_height * aspect_ratio, 1)
umap_width <- umap_panel_width + 46
umap_height <- umap_panel_height + 7

plot_umap <- function(field, colors, title) {
  p <- scop::CellDimPlot(
    object, reduction = "display", group.by = field,
    palcolor = colors, label = FALSE, seed = 11, show_stat = FALSE,
    raster = TRUE, raster.dpi = c(1600, 1600), pt.size = 2,
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    legend.title = "Cell type", theme_use = "theme_blank_axis",
    theme_args = list(text = element_text(family = "Arial")), combine = FALSE
  )[[1]]
  p <- order_dim_plot_cells(p, draw_order)
  stopifnot(identical(rownames(p$data), draw_order))
  p$layers <- Filter(function(layer) !inherits(layer$geom, "GeomCustomAnn"), p$layers)

  annotation_theme <- theme(
    text = element_text(family = "Arial", size = 6, face = "plain"),
    plot.title = element_text(size = 6.5, face = "plain", colour = "#1A1A1A"),
    plot.title.position = "plot",
    axis.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5.5),
    legend.key.height = grid::unit(2.2, "mm"),
    legend.key.width = grid::unit(2.2, "mm"),
    legend.spacing.x = grid::unit(0.8, "mm"),
    plot.margin = margin(1.5, 1.5, 3.5, 3.5, unit = "mm")
  )

  p + scale_colour_manual(name = "Cell type", values = colors, labels = identity) +
    theme_blank_axis(lab_size = 5.5, axis_lwd = 0.6,
                     xlab = "UMAP_1", ylab = "UMAP_2") +
    p$theme + coord_fixed(ratio = 1, clip = "off") +
    labs(title = title, subtitle = NULL) +
    guides(colour = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 1.3))) +
    annotation_theme
}

predicted_colors <- rcolors[names(rcolors) %in% unique(d$Predicted_CellType)]
ggsave(file.path(out, "fig4a.pdf"),
       plot_umap("Source_Label", qcolors, sprintf("Original labels (n = %d)", length(qcolors))),
       device = cairo_pdf, width = umap_width, height = umap_height, units = "mm")
ggsave(file.path(out, "fig4b.pdf"),
       plot_umap("Predicted_CellType", predicted_colors, sprintf("Transferred labels (n = %d)", length(predicted_colors))),
       device = cairo_pdf, width = umap_width, height = umap_height, units = "mm")

c <- fread(file.path(results_root, "full_source_prediction_counts.tsv"))
setnames(c, c("Source_CellClass", "Working_CellType", "Row_Percent"),
         c("Source_Label", "Predicted_CellType", "Source_Fraction"))
c[Source_Label %in% names(rename), Source_Label := unname(rename[Source_Label])]
c[, Source_Fraction := Source_Fraction / 100]

mat <- matrix(0, length(qcolors), length(unique(d$Predicted_CellType)),
              dimnames = list(names(qcolors), names(rcolors)[names(rcolors) %in% unique(d$Predicted_CellType)]))
mat[cbind(match(c$Source_Label, rownames(mat)), match(c$Predicted_CellType, colnames(mat)))] <- c$Source_Fraction
stopifnot(all(abs(rowSums(mat) - 1) < 1e-8), sum(c$Cells) == 1661798L)

# Label co-occurrence heatmap with compact layout, preserving % labels and clean grid lines.
cell_mm <- 3.2
ht <- Heatmap(
  mat,
  name = "Within source (%)",
  col = circlize::colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#6BAED6", "#08306B")),
  width = grid::unit(ncol(mat) * cell_mm, "mm"),
  height = grid::unit(nrow(mat) * cell_mm, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_rot = 50,
  row_names_gp = grid::gpar(fontsize = 5.2, fontfamily = "Arial"),
  column_names_gp = grid::gpar(fontsize = 5.2, fontfamily = "Arial"),
  row_names_max_width = grid::unit(38, "mm"),
  column_names_max_height = grid::unit(28, "mm"),
  row_labels = rownames(mat),
  left_annotation = rowAnnotation(
    simple_anno_size = grid::unit(2.0, "mm"),
    Source = rownames(mat),
    col = list(Source = qcolors),
    show_legend = FALSE,
    show_annotation_name = FALSE
  ),
  top_annotation = HeatmapAnnotation(
    simple_anno_size = grid::unit(2.0, "mm"),
    Predicted = colnames(mat),
    col = list(Predicted = rcolors),
    show_legend = FALSE,
    show_annotation_name = FALSE
  ),
  border = TRUE,
  rect_gp = grid::gpar(col = "#EDEDED", lwd = 0.3),
  heatmap_legend_param = list(
    at = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "25", "50", "75", "100"),
    legend_height = grid::unit(18, "mm"),
    grid_width = grid::unit(2.0, "mm"),
    title_gp = grid::gpar(fontsize = 5.2, fontface = "plain", fontfamily = "Arial"),
    labels_gp = grid::gpar(fontsize = 5.0, fontfamily = "Arial")
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    val <- mat[i, j]
    if (val >= 0.05) {
      grid::grid.text(
        sprintf("%.0f%%", 100 * val),
        x, y,
        gp = grid::gpar(
          fontsize = 4.2,
          fontfamily = "Arial",
          fontface = ifelse(val > 0.5, "bold", "plain"),
          col = ifelse(val > 0.5, "white", "#1A1A1A")
        )
      )
    }
  }
)

cairo_pdf(file.path(out, "fig4c.pdf"), width = 96 / 25.4, height = 70 / 25.4, family = "Arial")
draw(ht, heatmap_legend_side = "right", padding = grid::unit(c(1, 1, 1, 1), "mm"))
dev.off()

# Comparable donor concordance in a compact, wide layout
a <- fread(file.path(results_root, "comparable_donor_source_agreement.tsv"))
setnames(a, c("Source_Comparable", "Agreement"), c("Source_Label", "Concordance"))
short <- c("Neuron" = "Neurons", "Immune" = "Immune cells",
           "Vascular" = "Vascular cells", "Fibroblast" = "Perivascular fibroblasts")
a[Source_Label %in% names(short), Source_Label := unname(short[Source_Label])]

summary_stats <- a[, .(
  Mean = mean(Concordance),
  Median = median(Concordance),
  SD = sd(Concordance),
  N_Donors = .N,
  TotalCells = sum(Cells),
  CI_Lower = max(0, mean(Concordance) - qt(0.975, .N - 1) * sd(Concordance) / sqrt(.N)),
  CI_Upper = min(1, mean(Concordance) + qt(0.975, .N - 1) * sd(Concordance) / sqrt(.N))
), by = Source_Label]

summary_stats <- summary_stats[order(Mean)]
order_levels <- summary_stats$Source_Label
summary_stats[, Source_Label := factor(Source_Label, levels = order_levels)]
a[, Source_Label := factor(Source_Label, levels = order_levels)]

format_cells <- function(n) {
  if (n >= 1e5) sprintf("%.1fk", n / 1e3)
  else if (n >= 1e3) sprintf("%.1fk", n / 1e3)
  else as.character(n)
}
summary_stats[, CellStr := vapply(TotalCells, format_cells, character(1))]
summary_stats[, YLabel := sprintf("%s (%s)", Source_Label, CellStr)]
y_label_map <- setNames(summary_stats$YLabel, summary_stats$Source_Label)
summary_stats[, RightLabel := sprintf("%.1f%% (n = %d)", 100 * Mean, N_Donors)]

p_donor <- ggplot(summary_stats, aes(y = Source_Label)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = as.numeric(Source_Label) - 0.44,
                ymax = as.numeric(Source_Label) + 0.44,
                fill = Source_Label),
            inherit.aes = FALSE, alpha = 0.14, colour = NA) +
  geom_vline(xintercept = c(0, 0.5, 1.0),
             colour = "#E0E0E0", linewidth = 0.3) +
  geom_point(data = a, aes(x = Concordance, y = Source_Label),
             shape = 21, size = 1.3, stroke = 0.25,
             fill = "white", colour = "#4A4A4A", alpha = 0.85,
             position = position_jitter(width = 0, height = 0.16, seed = 2026)) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper, colour = Source_Label),
                orientation = "y", width = 0.28, linewidth = 0.45) +
  geom_point(aes(x = Mean, fill = Source_Label),
             shape = 23, size = 2.3, colour = "#1A1A1A", stroke = 0.4) +
  geom_text(aes(x = 1.04, label = RightLabel),
            hjust = 0, size = 2.05, family = "Arial", colour = "#2A2A2A") +
  scale_fill_manual(values = qcolors, guide = "none") +
  scale_colour_manual(values = qcolors, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1.0),
                     labels = c("0%", "50%", "100%"),
                     limits = c(-0.03, 1.38),
                     expand = c(0, 0)) +
  scale_y_discrete(labels = y_label_map, expand = expansion(add = 0.55)) +
  labs(x = "Donor concordance", y = NULL,
       title = "Donor concordance (605,701 cells from 26 donors)",
       subtitle = "Donor-equal mean (diamond) ± 95% CI") +
  theme_classic(base_family = "Arial", base_size = 6) +
  theme(
    plot.title = element_text(size = 6.2, face = "plain", colour = "#1A1A1A", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 5.0, colour = "#555555", margin = margin(b = 2)),
    axis.text.y = element_text(size = 5.3, colour = "#222222"),
    axis.text.x = element_text(size = 5.4, colour = "#222222"),
    axis.title.x = element_text(size = 5.8, colour = "#222222"),
    axis.line = element_line(linewidth = 0.3, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.3, colour = "#333333"),
    plot.margin = margin(1.5, 2, 1.5, 2, unit = "mm")
  )
ggsave(file.path(out, "fig4d.pdf"), p_donor, device = cairo_pdf, width = 85, height = 48, units = "mm")

# Assembly into standard manuscript double-column width (183 mm)
measure_panel_size <- function(path) {
  work <- tempfile("measure-")
  dir.create(work)
  cropped <- file.path(work, "crop.pdf")
  system2("pdfcrop", c("--margins", shQuote("2 2 2 2"), shQuote(path), shQuote(cropped)), stdout = FALSE)
  info <- system2("pdfinfo", shQuote(cropped), stdout = TRUE)
  size <- strsplit(sub("^Page size:\\s*", "", grep("^Page size:", info, value = TRUE)), "\\s+")[[1]]
  unlink(work, recursive = TRUE)
  c(width = as.numeric(size[1]) / 72 * 25.4, height = as.numeric(size[3]) / 72 * 25.4)
}

src <- normalizePath(file.path(out, "fig4a.pdf"))
pred <- normalizePath(file.path(out, "fig4b.pdf"))
lab <- normalizePath(file.path(out, "fig4c.pdf"))
don <- normalizePath(file.path(out, "fig4d.pdf"))

dim_a <- measure_panel_size(src)
dim_b <- measure_panel_size(pred)
dim_c <- measure_panel_size(lab)
dim_d <- measure_panel_size(don)

page_width <- 183
gap_x <- 6
gap_y <- 5

top_w <- (page_width - 8 - gap_x) / 2
top_h <- top_w * (dim_a[2] / dim_a[1])

bottom_w_c <- 87
bottom_h_c <- bottom_w_c * (dim_c[2] / dim_c[1])
bottom_w_d <- page_width - 8 - gap_x - bottom_w_c
bottom_h_d <- bottom_w_d * (dim_d[2] / dim_d[1])

bottom_row_h <- max(bottom_h_c, bottom_h_d)
page_height <- round(top_h + gap_y + bottom_row_h + 12, 1)

top_y <- page_height - 5
bottom_y <- top_y - top_h - gap_y

panels <- list(
  pdf_panel("A", 4, top_y, top_w, height = top_h, source = src, label_x = 1, label_y = top_y + 2),
  pdf_panel("B", 4 + top_w + gap_x, top_y, top_w, height = top_h, source = pred, label_x = 1 + top_w + gap_x, label_y = top_y + 2),
  pdf_panel("C", 4, bottom_y, bottom_w_c, height = bottom_h_c, source = lab, label_x = 1, label_y = bottom_y + 2),
  pdf_panel("D", 4 + bottom_w_c + gap_x, bottom_y, bottom_w_d, height = bottom_h_d, source = don, label_x = 1 + bottom_w_c + gap_x, label_y = bottom_y + 2)
)

assemble_pdf_figure(
  "fig4",
  c(page_width, page_height),
  panels,
  normalizePath(".")
)
