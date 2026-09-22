#!/usr/bin/env Rscript
# Figure 3: cluster markers, cell-type marker summary, source concordance.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
source("plotting/config.R")
studies <- fread(file.path(doc,"tables/source_concordance_dataset_harmonizable.tsv"))[is.finite(Agreement)]
d <- studies[, .(Dataset_Equal_Mean_Concordance=mean(Agreement),
 CI95_Lower=if(.N>1) max(0,mean(Agreement)-qt(.975,.N-1)*sd(Agreement)/sqrt(.N)) else NA_real_,
 CI95_Upper=if(.N>1) min(1,mean(Agreement)+qt(.975,.N-1)*sd(Agreement)/sqrt(.N)) else NA_real_),by=.(Main_CellType=Working_CellType)]
fwrite(d,"results/annotation/main_source_dataset_equal_summary.tsv",sep="\t")
d <- d[is.finite(d$Dataset_Equal_Mean_Concordance), ]
d <- d[order(d$Dataset_Equal_Mean_Concordance), ]
d$Main_CellType <- factor(d$Main_CellType, levels = rev(d$Main_CellType))
p <- ggplot(d, aes(Dataset_Equal_Mean_Concordance, Main_CellType)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
    ymin = as.numeric(Main_CellType) - .5, ymax = as.numeric(Main_CellType) + .5,
    fill = Main_CellType), inherit.aes = FALSE, alpha = .10, colour = NA) +
  scale_fill_manual(values = brainomics_celltype_colors, guide = "none") +
  geom_errorbar(aes(xmin = CI95_Lower, xmax = CI95_Upper), orientation = "y",
    width = .2, linewidth = .3, na.rm = TRUE) +
  geom_point(data=studies, aes(x=Agreement,y=Working_CellType),inherit.aes=FALSE,
    size=.6,shape=1,alpha=.6,colour="#777777") +
  geom_point(size = 1.2, colour = "#2E6E9E") +
  scale_x_continuous(limits = c(0, 1.03), breaks = c(0, .5, 1)) +
  labs(x = "Source-label concordance", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 6) +
  theme(axis.text.y = element_text(size = 5.5),
    plot.margin = margin(2, 2, 2, 2, "mm"))
ggsave("figures/fig3_concordance_layout_preview.pdf", p, device = grDevices::cairo_pdf,
  width = 73, height = 44, units = "mm", family = "Arial")
paths <- file.path("figures", c("cluster_marker_evidence_twocol_preview.pdf",
  "group_heatmap_flip.pdf", "fig3_concordance_layout_preview.pdf"))
widths <- c(176, 98, 73)
work <- tempfile("fig3-preview-"); dir.create(work)
heights <- vapply(seq_along(paths), function(i) {
  crop <- file.path(work, paste0(i, ".pdf"))
  system2("pdfcrop", c("--margins", shQuote("3 3 3 3"), shQuote(paths[i]), shQuote(crop)), stdout = FALSE)
  info <- system2("pdfinfo", shQuote(crop), stdout = TRUE)
  size <- strsplit(sub("^Page size:\\s*", "", grep("^Page size:", info, value = TRUE)), "\\s+")[[1]]
  widths[i] * as.numeric(size[3])/as.numeric(size[1])
}, numeric(1))
unlink(work, recursive = TRUE)
height <- heights[1] + max(heights[2:3]) + 15
tops <- c(height - 5, height - 10 - heights[1], height - 10 - heights[1])
xpos <- c(7, 7, 110)
assemble_pdf_figure("fig3", c(187, height),
  lapply(seq_along(paths), function(i) pdf_panel(LETTERS[i], xpos[i], tops[i], widths[i],
    source = normalizePath(paths[i]), label_x = xpos[i] - 5,
    label_y = tops[i] + 2)), normalizePath("."))
