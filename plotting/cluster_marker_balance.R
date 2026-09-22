#!/usr/bin/env Rscript
# Exploratory display of selected marker-set detection, not annotation confidence.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source("plotting/config.R")
setDTthreads(2)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
mapping <- a[,.(Cluster,CellType=Working_CellType,Cells)]
spec <- fread("results/annotation/display_markers.tsv")
evidence <- fread("results/annotation/full_cluster_marker_evidence.tsv")
stopifnot(!anyDuplicated(evidence[, .(Cluster, Gene)]))
# Controls are contextual competing-lineage/state evidence, not universal negatives.
sets <- fread("results/annotation/marker_role_sets.tsv")
sets <- rbind(sets[Direction == "Competing/context controls"], spec[,.(CellType=Marker_Group, Direction="Supporting markers", Gene)])
sets <- unique(sets, by=c("CellType","Gene"))
sets <- sets[!Gene %in% c("SLC1A2","SLC1A3","HDC","LUM","SLC32A1","TMEM119")]
negative <- split(sets[Direction == "Competing/context controls", Gene],
                  sets[Direction == "Competing/context controls", CellType])
stopifnot(!anyDuplicated(sets[, .(CellType, Gene)]), uniqueN(sets$Gene) >= 33L)
full_evidence <- evidence[Gene %in% unique(sets$Gene)]
stopifnot(nrow(full_evidence) == uniqueN(sets$Gene) * 75L)
fwrite(full_evidence, "results/annotation/full_validation_47_by_75.tsv", sep = "\t")
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
rows <- list(); headers <- list(); bands <- list(); cursor <- 0
for (ct in levels_type) {
  cl <- mapping[CellType == ct, Cluster]; cl <- cl[order(as.integer(sub("^C", "", cl)))]
  ng <- negative[[ct]]
  weighted_centers <- cumsum(nchar(ng) + 1.1) - nchar(ng) / 2
  negative_x <- if(length(ng)>1) 1 + (weighted_centers - min(weighted_centers)) / diff(range(weighted_centers)) * 10.6 else 6.3
  pg <- spec[Marker_Group == ct, Gene]
  support_spacing <- if (length(pg) >= 5L) 1.28 else if (length(pg) >= 4L) 1.2 else 1.65
  headers[[ct]] <- rbindlist(list(
    data.table(CellType = ct, Gene = ng, x = negative_x, y = cursor + 1.5),
    data.table(CellType = ct, Gene = pg, x = 15.35 + (seq_along(pg) - (length(pg) + 1) / 2) * support_spacing, y = cursor + 1.5)))
  rows[[ct]] <- data.table(CellType = ct, Cluster = cl,
    y = cursor + 1.45 + seq_along(cl))
  bands[[ct]] <- data.table(CellType = ct, y = cursor + .55,
    low = cursor + 2.0, high = cursor + 1.45 + length(cl) + .4)
  cursor <- cursor + 1.45 + length(cl) + .4
}
rows <- rbindlist(rows); headers <- rbindlist(headers); bands <- rbindlist(bands)
measured <- merge(measured, rows, by = c("CellType", "Cluster"))
measured <- merge(measured, headers[, .(CellType, Gene, x)], by = c("CellType", "Gene"))
stopifnot(uniqueN(measured$Cluster) == 75L,
  !anyDuplicated(measured[, .(Cluster, Direction, Gene)]))
fwrite(measured[, .(Cluster, CellType, Cells, Direction, Gene, Available_Cells,
  MeanLog, PctPositive, GeneZ)],
  "results/annotation/cluster_marker_balance_preview_source.tsv", sep = "\t")
p <- ggplot() +
  geom_rect(data = bands, aes(xmin = .4, xmax = 18.6, ymin = low, ymax = high),
    fill = "grey98", colour = "grey65", linewidth = .25) +
  geom_segment(data = rows, aes(x = .4, xend = 18.6, y = y, yend = y),
    colour = "grey90", linewidth = .15) +
  geom_segment(data = bands, aes(x = 12.5, xend = 12.5, y = low, yend = high),
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
  scale_x_continuous(limits = c(-.95, 19.1), expand = c(0, 0)) +
  annotate("text", x = 6.3, y = -.25, label = "Competing / context markers",
    family = "Arial", size = 2.8) +
  annotate("text", x = 15.1, y = -.25, label = "Supporting markers",
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
ggsave("figures/cluster_marker_evidence_preview.pdf", p,
  device = grDevices::cairo_pdf, width = 180, height = 360, units = "mm",
  family = "Arial", bg = "white")
message("75-row gene-level marker preview completed")

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
ggsave("figures/cluster_marker_evidence_twocol_preview.pdf", two_columns,
  device = grDevices::cairo_pdf, width = 220, height = 150, units = "mm",
  family = "Arial", bg = "white")
