#!/usr/bin/env Rscript
# Compose existing vector panels with the shared pdfcrop/XeLaTeX helper and
# write the PDF, SVG and night SVG of Figure 1.
# Keep the original Figure 1 narrative; split integration from annotation.
source("functions/utils.R")
repo <- normalizePath(".")
# Figure 1 uses equal gaps measured from the cropped panel dimensions.
work <- tempfile("fig1-layout-")
dir.create(work)
paths <- normalizePath(file.path("figures", c("fig1a_composition.pdf",
  "fig1b_age_coverage.pdf", "fig1c_workflow.pdf")))
# Match the assembly helper's crop before calculating equal visible spacing.
measure_panel_heights <- function(paths, widths) vapply(seq_along(paths), function(i) {
  cropped <- file.path(work, paste0("crop", i, ".pdf"))
  system2("pdfcrop", c("--margins", shQuote("3 3 3 3"), shQuote(paths[i]), shQuote(cropped)), stdout = FALSE)
  info <- system2("pdfinfo", shQuote(cropped), stdout = TRUE)
  size <- strsplit(sub("^Page size:\\s*", "", grep("^Page size:", info, value = TRUE)), "\\s+")[[1]]
  widths[i] * as.numeric(size[3]) / as.numeric(size[1])
}, numeric(1))
heights <- measure_panel_heights(paths, rep(168, length(paths)))
gap <- 5
page_height <- sum(heights) + (length(paths) - 1) * gap + 10
tops <- page_height - 5 - c(0, cumsum(head(heights, -1) + gap))
assemble_pdf_figure("fig1", c(183, page_height),
  lapply(seq_along(paths), function(i) pdf_panel(LETTERS[i], 8, tops[i], 168, source = paths[i])), repo)
# Keep the vector SVG and the inverted night SVG next to the PDF.
export_figure_assets(file.path(repo, "figures", "fig1.pdf"))
fig2_sources <- normalizePath(file.path("figures", c(
  "integration_umap.pdf", "fig2_validation_row.pdf", "fig2_annotation_row.pdf")))
fig2_heights <- measure_panel_heights(fig2_sources, rep(168, 3))
fig2_height <- sum(fig2_heights) + 2 * gap + 10
fig2_tops <- fig2_height - 5 - c(0, cumsum(head(fig2_heights, -1) + gap))
assemble_pdf_figure("fig2", c(183, fig2_height), list(
  pdf_panel("A", 8, fig2_tops[1], 168, source = fig2_sources[1]),
  pdf_panel("BCDE", 8, fig2_tops[2], 168, source = fig2_sources[2], label = ""),
  pdf_panel("FG", 8, fig2_tops[3], 168, source = fig2_sources[3], label = "")
), repo)
unlink(work, recursive = TRUE)
# fig3 is assembled by plotting/fig3_annotation_panels.R and fig4 by the EGAD
# step; this script writes the Fig. 1 and Fig. 2 composites.
source("plotting/fig4.R")
