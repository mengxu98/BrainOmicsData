#!/usr/bin/env Rscript
# Exploratory display of selected marker-set detection, not annotation confidence.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source("functions/config.R")
setDTthreads(2)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
mapping <- a[,.(Cluster,CellType=Working_CellType,Cells)]
spec <- fread("results/annotation/display_markers.tsv")
# Measured genes omitted from the 33-gene display: role-set support plus
# canonical markers that were filtered out of the preview.
spec <- rbind(spec, data.table(
  Marker_Group = c(
    "Astrocytes", "Astrocytes", "Astrocytes",
    "Mural cells", "Fibroblasts", "Fibroblasts",
    "Excitatory neurons", "Excitatory neurons", "Excitatory neurons",
    "Inhibitory neurons", "Inhibitory neurons",
    "Microglia", "Microglia", "Neural progenitors", "Oligodendrocytes"),
  Gene = c(
    "GJA1", "SLC1A2", "SLC1A3",
    "PDGFRB", "COL1A2", "LUM",
    "SATB2", "TBR1", "STMN2",
    "SLC32A1", "STMN2",
    "CSF1R", "TMEM119", "FABP7", "MOBP")))
spec <- unique(spec, by = c("Marker_Group", "Gene"))
evidence <- fread("results/annotation/full_cluster_marker_evidence.tsv")
stopifnot(!anyDuplicated(evidence[, .(Cluster, Gene)]))
# Controls are contextual competing-lineage/state evidence, not universal negatives.
sets <- fread("results/annotation/marker_role_sets.tsv")
sets <- rbind(sets[Direction == "Competing/context controls"], spec[,.(CellType=Marker_Group, Direction="Supporting markers", Gene)])
sets <- unique(sets, by=c("CellType","Gene"))
# One displayed gene per competing lineage. The role table itself is unchanged.
sets <- sets[!data.table(
  CellType = c(rep("Astrocytes", 5), rep("Microglia", 5),
    rep("Neural progenitors", 2), rep("Oligodendrocyte progenitor cells", 3),
    rep("Differentiating oligodendrocytes", 3), rep("Oligodendrocytes", 3),
    rep("Mural cells", 2), rep("Fibroblasts", 2), rep("Lymphocytes", 2)),
  Gene = c("PECAM1", "CSF1R", "MBP", "TP73", "DNAH11",
    "MBP", "PLP1", "MOBP", "MRC1", "LYVE1",
    "ALDH1L1", "MBP",
    "CSF1R", "PECAM1", "GAD2",
    "CSF1R", "PECAM1", "GAD2",
    "CSF1R", "PECAM1", "GAD2",
    "VWF", "COL1A2",
    "VWF", "ACTA2",
    "CSF1R", "C1QA")
), on = c("CellType", "Gene")]
negative <- split(sets[Direction == "Competing/context controls", Gene],
                  sets[Direction == "Competing/context controls", CellType])
stopifnot(!anyDuplicated(sets[, .(CellType, Gene)]), uniqueN(sets$Gene) >= 33L)
full_evidence <- evidence[Gene %in% unique(sets$Gene)]
stopifnot(nrow(full_evidence) == uniqueN(sets$Gene) * 75L)
fwrite(full_evidence, "results/annotation/full_validation_marker_by_75.tsv", sep = "\t")
measured <- merge(merge(mapping, sets, by = "CellType", allow.cartesian = TRUE),
  evidence[, .(Cluster, Gene, Available_Cells, MeanLog, PctPositive)],
  by = c("Cluster", "Gene"), all.x = TRUE)
stopifnot(all(is.finite(measured$PctPositive)),
  all(measured$PctPositive >= 0 & measured$PctPositive <= 1),
  all(measured$Available_Cells > 0 & measured$Available_Cells <= measured$Cells))
# One dot per measured gene; no averaging across marker sets.
levels_type <- names(brainomics_celltype_colors)
levels_type <- levels_type[levels_type %in% mapping$CellType]
# Move the short lymphocyte block to the left to balance both columns.
levels_type <- c(levels_type[1:5], "Lymphocytes", setdiff(levels_type[6:12], "Lymphocytes"))
evidence[, GeneZ := as.numeric(scale(MeanLog)), by = Gene]
measured <- merge(measured, evidence[, .(Cluster, Gene, GeneZ)],
  by = c("Cluster", "Gene"), all.x = TRUE)
stopifnot(all(is.finite(measured$GeneZ)))
# One x unit is the inner panel millimetre at the saved column width.
inner_mm_per_x <- 5.39
label_gap_mm <- 2.0
label_device <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(label_device, width = 4, height = 2, family = "Arial")
on.exit(grDevices::dev.off(), add = TRUE)
grid::grid.newpage()
label_gp <- grid::gpar(fontfamily = "Arial", fontface = "italic", fontsize = 1.85 * 72.27 / 25.4)
gene_width_mm <- function(genes) vapply(genes, function(gene) {
  grid::convertWidth(grid::grobWidth(grid::textGrob(gene, gp = label_gp)), "mm", valueOnly = TRUE)
}, numeric(1))
rows <- list(); headers <- list(); bands <- list(); cursor <- 0
competing_left <- 0.85
competing_gap_mm <- 1.6
negative_x <- list()
competing_right <- competing_left
for (ct in levels_type) {
  ng <- negative[[ct]]
  widths <- gene_width_mm(ng)
  xs <- numeric(length(ng))
  xs[1] <- competing_left + widths[1] / 2 / inner_mm_per_x
  if (length(ng) > 1L) {
    for (i in 2:length(ng))
      xs[i] <- xs[i - 1] + (widths[i - 1] + widths[i]) / 2 / inner_mm_per_x + competing_gap_mm / inner_mm_per_x
  }
  negative_x[[ct]] <- xs
  competing_right <- max(competing_right, xs[length(ng)] + widths[length(ng)] / 2 / inner_mm_per_x)
}
support_left <- competing_right + 0.7
support_right <- support_left
divider_x <- competing_right + 0.35
for (ct in levels_type) {
  cl <- mapping[CellType == ct, Cluster]; cl <- cl[order(as.integer(sub("^C", "", cl)))]
  ng <- negative[[ct]]
  pg <- spec[Marker_Group == ct, Gene]
  widths <- gene_width_mm(pg)
  support_x <- numeric(length(pg))
  support_x[1] <- support_left + widths[1] / 2 / inner_mm_per_x
  if (length(pg) > 1L) {
    for (i in 2:length(pg))
      support_x[i] <- support_x[i - 1] + (widths[i - 1] + widths[i]) / 2 / inner_mm_per_x + label_gap_mm / inner_mm_per_x
  }
  support_right <- max(support_right, support_x[length(pg)] + widths[length(pg)] / 2 / inner_mm_per_x)
  headers[[ct]] <- rbindlist(list(
    data.table(CellType = ct, Gene = ng, x = negative_x[[ct]], y = cursor + 1.5),
    data.table(CellType = ct, Gene = pg, x = support_x, y = cursor + 1.5)))
  rows[[ct]] <- data.table(CellType = ct, Cluster = cl,
    y = cursor + 1.45 + seq_along(cl))
  bands[[ct]] <- data.table(CellType = ct, y = cursor + .55,
    low = cursor + 2.0, high = cursor + 1.45 + length(cl) + .4)
  cursor <- cursor + 1.45 + length(cl) + .4
}
rows <- rbindlist(rows); headers <- rbindlist(headers); bands <- rbindlist(bands)
panel_right <- support_right + 0.28
x_max <- panel_right + 0.45
support_headers <- headers[x >= support_left - 0.01]
support_headers[, width_mm := gene_width_mm(Gene)]
support_gaps <- support_headers[, .(
  gap_mm = (x - shift(x)) * inner_mm_per_x - (width_mm + shift(width_mm)) / 2
), by = CellType]
stopifnot(all(support_gaps[is.finite(gap_mm), gap_mm] >= label_gap_mm - 0.05))
measured <- merge(measured, rows, by = c("CellType", "Cluster"))
measured <- merge(measured, headers[, .(CellType, Gene, x)], by = c("CellType", "Gene"))
stopifnot(uniqueN(measured$Cluster) == 75L,
  !anyDuplicated(measured[, .(Cluster, Direction, Gene)]))
fwrite(measured[, .(Cluster, CellType, Cells, Direction, Gene, Available_Cells,
  MeanLog, PctPositive, GeneZ)],
  "results/annotation/cluster_marker_balance_preview_source.tsv", sep = "\t")
p <- ggplot() +
  geom_rect(data = bands, aes(xmin = .4, xmax = panel_right, ymin = low, ymax = high),
    fill = "grey98", colour = "grey65", linewidth = .25) +
  geom_segment(data = rows, aes(x = .4, xend = panel_right, y = y, yend = y),
    colour = "grey90", linewidth = .15) +
  geom_segment(data = bands, aes(y = low, yend = high), x = divider_x, xend = divider_x,
    colour = "grey65", linewidth = .25) +
  geom_point(data = measured, aes(x = x, y = y, size = PctPositive, fill = GeneZ),
    shape = 21, colour = "grey20", stroke = .15) +
  geom_text(data = rows, aes(x = .1, y = y, label = Cluster),
    hjust = 1, family = "Arial", size = 2.1) +
  geom_text(data = headers, aes(x = x, y = y, label = Gene),
    family = "Arial", fontface = "italic", size = 1.85) +
  geom_text(data = bands, aes(x = -.35, y = y, label = CellType, colour = CellType),
    hjust = 0, family = "Arial", size = 2.4) +
  scale_colour_manual(values = brainomics_celltype_colors, guide = "none") +
  scale_size_area(name = "Detected (%)", limits = c(0, 1), max_size = 2.0,
    breaks = c(.2, .5, 1), labels = c("20%", "50%", "100%")) +
  scale_fill_gradientn(name = "Z-score",
    colours = c("#3288BD", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D53E4F"),
    limits = c(-2, 2), oob = scales::squish, breaks = c(-2, 0, 2),
    labels = c("-2", "0", "2")) +
  scale_y_reverse(limits = c(-.75, cursor + .25), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-.95, x_max), expand = c(0, 0)) +
  annotate("text", x = (competing_left + competing_right) / 2, y = -.25, label = "Competing markers",
    family = "Arial", size = 2.8) +
  annotate("text", x = (support_left + support_right) / 2, y = -.25, label = "Supporting markers",
    family = "Arial", size = 2.8) +
  guides(size = guide_legend(order = 1, override.aes = list(fill = "grey55")),
    fill = guide_colourbar(order = 2, barwidth = grid::unit(16, "mm"),
      barheight = grid::unit(2, "mm"))) +
  theme_void(base_family = "Arial", base_size = 7) +
  theme(legend.position = "bottom", legend.box = "horizontal",
    plot.title = element_text(size = 10, face = "plain", margin = margin(b = 6)),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
    legend.key.width = grid::unit(3, "mm"), legend.key.height = grid::unit(3, "mm"),
    legend.spacing.x = grid::unit(1, "mm"),
    legend.box.spacing = grid::unit(1, "mm"),
    plot.margin = margin(3, 3, 3, 3, "mm"))

# Two-column layout, keeping each cell type intact; balance occupied row space.
split_at <- 6L
column_types <- list(levels_type[seq_len(split_at)], levels_type[-seq_len(split_at)])
column_offsets <- vapply(column_types, function(ct) min(bands[CellType %in% ct, y]) - .55, numeric(1))
column_height <- max(vapply(seq_along(column_types), function(i)
  max(bands[CellType %in% column_types[[i]], high]) - column_offsets[i], numeric(1)))
column_plots <- lapply(seq_along(column_types), function(i) {
  q <- unserialize(serialize(p, NULL))
  q$layers <- lapply(q$layers, function(layer) {
    d <- layer$data
    if (is.data.frame(d) && "CellType" %in% names(d)) {
      d <- as.data.table(copy(d))[CellType %in% column_types[[i]]]
      for (field in intersect(c("y", "low", "high"), names(d)))
        set(d, j = field, value = d[[field]] - column_offsets[i])
      # Keep original gene positions unchanged between columns.
      layer$data <- d
    }
    layer
  })
  q + scale_y_reverse(limits = c(-.75, column_height + .15), expand = c(0, 0)) +
    labs(title = NULL) + theme(plot.margin = margin(.5, 1, .5, 1, "mm"))
})
two_columns <- patchwork::wrap_plots(column_plots, nrow = 1, guides = "collect")
two_columns <- two_columns & theme(legend.position = "bottom")
page_width_mm <- 2 * ((x_max - (-.95)) * inner_mm_per_x + 2)
ggsave("figures/fig3a.pdf", two_columns,
  device = grDevices::cairo_pdf, width = page_width_mm, height = 150, units = "mm",
  family = "Arial", bg = "white")
message("75-row two-column marker figure completed")
