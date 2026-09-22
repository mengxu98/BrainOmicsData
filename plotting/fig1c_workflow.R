#!/usr/bin/env Rscript
# Figure 1C: the analysis workflow from source collection to cell-type annotation.
suppressPackageStartupMessages({library(ggplot2)})
source("functions/utils.R")

panel_dir <- "figures"
asset_dir <- "results/annotation"
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)
workflow <- data.frame(
  step = factor(seq_len(6L)),
  label = c(
    "Source\ncollection", "Metadata\ncuration", "Age/region\nmapping",
    "Common-count\nintegration", "Method\nevaluation",
    "Cell type\nannotation"
  ),
  detail = c(
    "22 datasets", "286 donors\n448 specimens\n3,744 libraries",
    "15 intervals\n48 regions", "2,602,031 cells\n50 dimensions",
    "Raw, scVI\nHarmony, RPCA", "75 clusters\n12 cell types"
  ),
  x = seq(1, 11, by = 2), y = 1, stringsAsFactors = FALSE
)
arrows <- data.frame(
  x = workflow$x[-nrow(workflow)], xend = workflow$x[-1L],
  y = 1, yend = 1
)
box_w <- 1.46
box_h <- 0.92 * 5 / 6
box_r <- 0.16 * 5 / 6
rounded_box <- function(xmin, xmax, ymin, ymax, r, n = 20L) {
  r <- min(r, (xmax - xmin) / 2, (ymax - ymin) / 2)
  arc <- function(cx, cy, a0, a1) {
    a <- seq(a0, a1, length.out = n)
    data.frame(px = cx + r * cos(a), py = cy + r * sin(a))
  }
  rbind(
    arc(xmax - r, ymax - r, 0, pi / 2),
    arc(xmin + r, ymax - r, pi / 2, pi),
    arc(xmin + r, ymin + r, pi, 3 * pi / 2),
    arc(xmax - r, ymin + r, 3 * pi / 2, 2 * pi)
  )
}
boxes <- do.call(rbind, lapply(seq_len(nrow(workflow)), function(i) {
  path <- rounded_box(
    workflow$x[[i]] - box_w / 2, workflow$x[[i]] + box_w / 2,
    workflow$y[[i]] - box_h / 2, workflow$y[[i]] + box_h / 2,
    box_r
  )
  path$step <- workflow$step[[i]]
  path
}))
fig1c <- ggplot(workflow, aes(x, y)) +
  geom_segment(
    data = arrows,
    aes(x = x + box_w / 2, xend = xend - box_w / 2, y = y, yend = yend),
    linewidth = 0.22, colour = "#56616D",
    arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")
  ) +
  geom_polygon(
    data = boxes,
    aes(px, py, group = step),
    fill = "#F2F2F2", colour = "#A4B4BF", linewidth = 0.2, inherit.aes = FALSE
  ) +
  geom_text(aes(label = label),
    nudge_y = 0.20 * 5 / 6, size = 2.15,
    fontface = "plain", lineheight = 0.90, colour = "#111827"
  ) +
  geom_text(aes(label = detail),
    nudge_y = -0.21 * 5 / 6, size = 1.72,
    lineheight = 0.90, colour = "#52616B"
  ) +
  coord_cartesian(
    xlim = c(0.18, 11.82),
    ylim = c(1 - box_h / 2 - 0.10, 1 + box_h / 2 + 0.10),
    clip = "on"
  ) +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(1.2, 1.5, 1.2, 1.5))

write_tsv(workflow, file.path(asset_dir, "figure1_workflow_steps.tsv"))

ggplot2::ggsave(
  file.path(panel_dir, "fig1c_workflow.pdf"), fig1c,
  device = grDevices::cairo_pdf, width = 167, height = 20, units = "mm",
  family = "Arial", bg = "white"
)
