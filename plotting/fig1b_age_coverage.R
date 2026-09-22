#!/usr/bin/env Rscript
# Figure 1B: reported-age coverage with stage totals and intervals.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("functions/dataset_metadata.R")
analysis_dir <- Sys.getenv(
  "BRAINOMICS_ANALYSIS_DIR",
  unset = file.path(
    "outputs", "reprocessing_plan_20260913", "execution_20260914_34387",
    "analysis_from_scratch_20260917"
  )
)
summary_dir <- Sys.getenv(
  "BRAINOMICS_REFERENCE_SUMMARY",
  unset = file.path(analysis_dir, "07_downstream", "revision_20260918", "reference_summary")
)
ages <- fread(file.path(summary_dir, "reported_age_summary.tsv"))
specimens <- fread(file.path(summary_dir, "reported_age_sex_specimen_summary.tsv"))
stages <- fread(file.path(summary_dir, "age_interval_summary.tsv"))
setorder(ages, Age_Display_Order)
ages[, x := as.numeric(seq_len(.N))]
ages[, AgeIntervalID := dataset_age_interval(
  Age_num, ifelse(Unit == "Years", "years", Unit)
)]
specimens[, x := ages$x[match(Age_Label, ages$Age_Label)]]
specimens[, AgeIntervalID := ages$AgeIntervalID[match(Age_Label, ages$Age_Label)]]
bands <- ages[, .(left = min(x) - 0.5, right = max(x) + 0.5,
  Cells = sum(Cells)), by = AgeIntervalID]
bands[, mid := (left + right) / 2]
specimen_totals <- specimens[, .(Specimens = sum(Specimens)), by = AgeIntervalID]
bands[, Specimens := specimen_totals$Specimens[match(AgeIntervalID, specimen_totals$AgeIntervalID)]]
bands[, AgeRange := stages$AgeRange[match(AgeIntervalID, stages$AgeIntervalID)]]
bands[, Interval_Label := sub(" (PCW|years)$", "", AgeRange)]
bands[, Interval_Label := sub("+Inf", "∞", Interval_Label, fixed = TRUE)]
# The existing export verifies one age and sex per known specimen. Summing
# these disjoint age-sex strata therefore counts specimens, not donors.
stopifnot(!anyNA(ages$AgeIntervalID), !anyNA(specimens$x),
  identical(bands$Cells, stages$Cells[match(bands$AgeIntervalID, stages$AgeIntervalID)]),
  sum(bands$Specimens) == 448L)
# Separate prenatal and postnatal age categories without hiding observations.
age_gap <- 0.6
postnatal_stages <- paste0("S", 8:15)
ages[AgeIntervalID %in% postnatal_stages, x := x + age_gap]
specimens[AgeIntervalID %in% postnatal_stages, x := x + age_gap]
bands[AgeIntervalID %in% postnatal_stages,
  `:=`(left = left + age_gap, right = right + age_gap, mid = mid + age_gap)]
x_span <- nrow(ages) + age_gap
stage_colors <- setNames(c(
  colorRampPalette(c("#0AA344", "#006D87"))(7),
  colorRampPalette(c("#2B73AF", "#003D74"))(8)
), paste0("S", 1:15))
sex_colors <- c(Female = "#B24AA6", Male = "#28327F", "Not reported" = "#8A8A8A")
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
base <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  axis.text = element_text(colour = "#303030"),
  panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.2),
  axis.line = element_line(linewidth = 0.3),
  plot.margin = margin(1, 3, 1, 3)
)
x_scale <- scale_x_continuous(limits = c(0.5, x_span + 0.5), expand = c(0, 0))
background <- geom_rect(data = bands, aes(xmin = left, xmax = right,
  ymin = -Inf, ymax = Inf, fill = AgeIntervalID), inherit.aes = FALSE,
  alpha = 0.12)
fill_scale <- scale_fill_manual(values = stage_colors, guide = "none")
boundaries <- geom_vline(xintercept = head(bands$right, -1), colour = "white", linewidth = 0.25)
strip_theme <- theme_void(base_family = "Arial") + theme(plot.margin = margin(1, 3, 1, 3))
unit_groups <- data.table(
  left = c(bands$left[1], bands$left[8]),
  right = c(bands$right[7], bands$right[15]),
  label = c("Prenatal (PCW)", "Postnatal (years)")
)
unit_groups[, mid := (left + right) / 2]
top <- ggplot(unit_groups) +
  geom_rect(data = bands, aes(xmin = left, xmax = right,
    ymin = 0, ymax = 1, fill = AgeIntervalID)) +
  geom_text(aes(mid, 0.5, label = label), colour = "white", size = 2.2) +
  fill_scale +
  x_scale + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + strip_theme +
  theme(plot.margin = margin(1, 3, 0, 3))
cells <- ggplot(ages, aes(x, Cells)) + background + fill_scale + boundaries +
  geom_col(width = 0.72, fill = "black") + x_scale +
  scale_y_continuous(labels = fmt, limits = c(0, max(ages$Cells) * 1.25),
    breaks = c(0, 50000, 100000), expand = c(0, 0)) +
  labs(x = NULL, y = "Cells/nuclei") + base +
  geom_segment(data = unit_groups, aes(x = left, xend = right, y = 0, yend = 0),
    inherit.aes = FALSE, linewidth = 0.3) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), plot.margin = margin(0, 3, 1, 3))
sex_totals <- specimens[, .(Specimens = sum(Specimens)), by = Reported_Sex]
sex_labels <- setNames(paste0(sex_totals$Reported_Sex, " (n = ", sex_totals$Specimens, ")"), sex_totals$Reported_Sex)
points <- ggplot(specimens, aes(x, Specimens, size = Specimens, colour = Reported_Sex)) +
  background + fill_scale + boundaries + geom_point(position = position_dodge(width = 0.7), alpha = 0.9) +
  geom_point(aes(group = Reported_Sex), position = position_dodge(width = 0.7),
    shape = 21, fill = NA, colour = "white", alpha = 0.55, stroke = 0.18,
    show.legend = FALSE) +
  geom_text(data = bands, aes(mid, max(specimens$Specimens) + 3.5, label = fmt(Specimens)),
    inherit.aes = FALSE, size = 2.15, family = "Arial") +
  scale_x_continuous(limits = c(0.5, x_span + 0.5), expand = c(0, 0),
    breaks = NULL) +
  scale_y_reverse(limits = c(max(specimens$Specimens) + 5.5, 0),
    breaks = c(5, 10, 15, 20), expand = c(0, 0)) +
  scale_colour_manual(values = sex_colors, labels = sex_labels, name = "Reported sex") +
  scale_size_continuous(range = c(0.7, 2.6), guide = "none") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2))) +
  labs(x = NULL, y = "Specimens") + base +
  geom_segment(data = unit_groups,
    aes(x = left, xend = right, y = max(specimens$Specimens) + 5.5,
      yend = max(specimens$Specimens) + 5.5),
    inherit.aes = FALSE, linewidth = 0.3) +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), plot.margin = margin(1, 3, 0.425, 3),
    legend.position = "bottom", legend.text = element_text(size = 6),
    legend.title = element_text(size = 6), legend.key.size = grid::unit(3, "mm"))
bottom <- ggplot(bands) +
  geom_rect(aes(xmin = left + 0.25, xmax = right - 0.25,
    ymin = 0.83, ymax = 1.15, fill = AgeIntervalID)) +
  geom_text(aes(mid, 0.99, label = AgeIntervalID,
    angle = ifelse(right - left <= 3 & nchar(AgeIntervalID) > 2, 90, 0)),
    colour = "white", size = 1.9) +
  geom_text(aes(mid, 0.76, label = Interval_Label), angle = 90, hjust = 1, size = 2.0) +
  fill_scale + x_scale +
  scale_y_continuous(limits = c(0, 1.15), expand = c(0, 0)) + strip_theme +
  theme(plot.margin = margin(0, 3, 1, 3))
# Measure text and available panel width on the export device. Rotation is
# determined by each stage's physical width, not a hard-coded list of stages.
figure_width_mm <- 183
figure_height_mm <- 82
measurement_file <- tempfile(fileext = ".pdf")
cairo_pdf(measurement_file, width = figure_width_mm / 25.4, height = figure_height_mm / 25.4, family = "Arial")
body_grob <- patchwork::patchworkGrob(cells / points)
panel_width_mm <- figure_width_mm - grid::convertWidth(sum(body_grob$widths), "mm", valueOnly = TRUE)
text_width_mm <- vapply(fmt(bands$Cells), function(label) {
  grid::convertWidth(grid::grobWidth(grid::textGrob(label,
    gp = grid::gpar(fontfamily = "Arial", fontsize = 2.1 * 72.27 / 25.4))), "mm", valueOnly = TRUE)
}, numeric(1))
invisible(dev.off())
unlink(measurement_file)
bands[, count_angle := ifelse(text_width_mm + 1 > (right - left) / x_span * panel_width_mm, 90, 0)]
cells <- cells + geom_text(data = bands,
  aes(mid, max(ages$Cells) * 1.18, label = fmt(Cells), angle = count_angle,
    hjust = ifelse(count_angle == 90, 1, 0.5), vjust = ifelse(count_angle == 90, 0.5, 1)),
  inherit.aes = FALSE, family = "Arial", size = 2.1)
preview <- (top / cells / points / bottom) +
  plot_layout(heights = c(0.25, 2.5, 1.8, 0.9), guides = "collect") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(0, "mm"))
output <- "figures/fig1b_age_coverage.pdf"
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
ggsave(output, preview, width = figure_width_mm, height = figure_height_mm, units = "mm", device = cairo_pdf, family = "Arial")
message("Figure 1B: ", output)
