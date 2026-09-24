#!/usr/bin/env Rscript
# Figure 3: cluster markers, cell-type marker summary, source concordance.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
source("functions/config.R")
studies <- fread(file.path(doc,"tables/source_concordance_dataset_harmonizable.tsv"))[is.finite(Agreement)]
d <- studies[, .(Dataset_Equal_Mean_Concordance=mean(Agreement),
 CI95_Lower=if(.N>1) max(0,mean(Agreement)-qt(.975,.N-1)*sd(Agreement)/sqrt(.N)) else NA_real_,
 CI95_Upper=if(.N>1) min(1,mean(Agreement)+qt(.975,.N-1)*sd(Agreement)/sqrt(.N)) else NA_real_),by=.(Main_CellType=Working_CellType)]
fwrite(d,"results/annotation/main_source_dataset_equal_summary.tsv",sep="\t")
d <- d[is.finite(d$Dataset_Equal_Mean_Concordance), ]
d <- d[order(d$Dataset_Equal_Mean_Concordance), ]
d$Main_CellType <- factor(d$Main_CellType, levels = rev(d$Main_CellType))
studies$Working_CellType <- factor(studies$Working_CellType, levels = levels(d$Main_CellType))
p <- ggplot(d, aes(Dataset_Equal_Mean_Concordance, Main_CellType)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
    ymin = as.numeric(Main_CellType) - .42, ymax = as.numeric(Main_CellType) + .42,
    fill = Main_CellType), inherit.aes = FALSE, alpha = .14, colour = NA) +
  geom_errorbar(aes(xmin = CI95_Lower, xmax = CI95_Upper, colour = Main_CellType),
    orientation = "y", width = .28, linewidth = .45, na.rm = TRUE) +
  geom_point(data = studies, aes(x = Agreement, y = Working_CellType),
    inherit.aes = FALSE, shape = 21, size = 1.15, stroke = .25,
    fill = "white", colour = "#4A4A4A",
    position = position_jitter(width = 0, height = .16, seed = 1)) +
  geom_point(aes(fill = Main_CellType), shape = 21, size = 1.9,
    colour = "black", stroke = .3) +
  scale_fill_manual(values = brainomics_celltype_colors, guide = "none") +
  scale_colour_manual(values = brainomics_celltype_colors, guide = "none") +
  scale_x_continuous(limits = c(-.03, 1.04), breaks = c(0, .5, 1),
    expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = .55)) +
  labs(x = "Source-label concordance", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 6) +
  theme(axis.text.y = element_text(size = 5.4, colour = "#222222"),
    axis.text.x = element_text(size = 5.4, colour = "#222222"),
    axis.title.x = element_text(size = 6, colour = "#222222"),
    axis.line = element_line(linewidth = .25, colour = "black"),
    axis.ticks = element_line(linewidth = .25, colour = "black"),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))
ggsave("figures/fig3c.pdf", p, device = grDevices::cairo_pdf,
  width = 78, height = 52, units = "mm", family = "Arial")
paths <- file.path("figures", c("fig3a.pdf",
  "fig3b.pdf", "fig3c.pdf"))
work <- tempfile("fig3-layout-"); dir.create(work)
aspect <- vapply(paths, function(path) {
  crop <- file.path(work, basename(path))
  system2("pdfcrop", c("--margins", shQuote("3 3 3 3"), shQuote(path), shQuote(crop)), stdout = FALSE)
  info <- system2("pdfinfo", shQuote(crop), stdout = TRUE)
  size <- strsplit(sub("^Page size:\\s*", "", grep("^Page size:", info, value = TRUE)), "\\s+")[[1]]
  as.numeric(size[3]) / as.numeric(size[1])
}, numeric(1))
unlink(work, recursive = TRUE)
# B and C share one row height; the flatter heatmap takes the wider share.
content_width <- 176
bottom_gap <- 5
width_c <- (content_width - bottom_gap) / (aspect[3] / aspect[2] + 1)
width_b <- content_width - bottom_gap - width_c
widths <- c(content_width, width_b, width_c)
heights <- widths * aspect
height <- heights[1] + heights[2] + 15
tops <- c(height - 5, height - 10 - heights[1], height - 10 - heights[1])
xpos <- c(7, 7, 7 + width_b + bottom_gap)
assemble_pdf_figure("fig3", c(187, height),
  lapply(seq_along(paths), function(i) pdf_panel(LETTERS[i], xpos[i], tops[i], widths[i],
    source = normalizePath(paths[i]), label_x = xpos[i] - 5,
    label_y = tops[i] + 2)), normalizePath("."))
