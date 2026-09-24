#!/usr/bin/env Rscript
# Compose existing vector panels with the shared pdfcrop/XeLaTeX helper and
# write the PDF, SVG and night SVG of Figure 1.
# Keep the original Figure 1 narrative; split integration from annotation.
source("functions/utils.R")
repo <- normalizePath(".")
target <- Sys.getenv("BRAINOMICS_ASSEMBLE", unset = "all")
stopifnot(target %in% c("fig1", "fig2", "all"))
# Figure 1 uses equal gaps measured from the cropped panel dimensions.
work <- tempfile("fig1-layout-")
dir.create(work)
gap <- 5
# Match the assembly helper's crop before calculating equal visible spacing.
measure_panel_heights <- function(paths, widths) vapply(seq_along(paths), function(i) {
  cropped <- file.path(work, paste0("crop", i, ".pdf"))
  system2("pdfcrop", c("--margins", shQuote("3 3 3 3"), shQuote(paths[i]), shQuote(cropped)), stdout = FALSE)
  info <- system2("pdfinfo", shQuote(cropped), stdout = TRUE)
  size <- strsplit(sub("^Page size:\\s*", "", grep("^Page size:", info, value = TRUE)), "\\s+")[[1]]
  widths[i] * as.numeric(size[3]) / as.numeric(size[1])
}, numeric(1))
if (target %in% c("fig1", "all")) {
paths <- normalizePath(file.path("figures", c("fig1a.pdf",
  "fig1b.pdf", "fig1c.pdf")))
heights <- measure_panel_heights(paths, rep(168, length(paths)))
page_height <- sum(heights) + (length(paths) - 1) * gap + 10
tops <- page_height - 5 - c(0, cumsum(head(heights, -1) + gap))
assemble_pdf_figure("fig1", c(183, page_height),
  lapply(seq_along(paths), function(i) pdf_panel(LETTERS[i], 8, tops[i], 168, source = paths[i])), repo)
# Keep the vector SVG and the inverted night SVG next to the PDF.
export_figure_assets(file.path(repo, "figures", "fig1.pdf"))
}
if (target %in% c("fig2", "all")) {
fig2_sources <- normalizePath(file.path("figures", c(
  "fig2a.pdf", "fig2b.pdf",
  "fig2c.pdf", "fig2d.pdf",
  "fig2e.pdf", "fig2f.pdf")))
fig2_widths <- c(176, rep(41, 4), 176)
fig2_heights <- measure_panel_heights(fig2_sources, fig2_widths)
fig2_row_heights <- c(fig2_heights[1], max(fig2_heights[2:5]),
                      fig2_heights[6])
fig2_height <- sum(fig2_row_heights) + 2 * gap + 10
fig2_tops <- fig2_height - 5 - c(0, cumsum(head(fig2_row_heights, -1) + gap))
fig2_x <- c(7, 7, 52, 97, 142, 7)
fig2_y <- fig2_tops[c(1, 2, 2, 2, 2, 3)]
assemble_pdf_figure("fig2", c(187, fig2_height),
  lapply(seq_along(fig2_sources), function(i) pdf_panel(
    LETTERS[i], fig2_x[i], fig2_y[i], fig2_widths[i],
    source = fig2_sources[i], label_x = fig2_x[i] - 5,
    label_y = fig2_y[i] + 2)), repo)
}
unlink(work, recursive = TRUE)
# fig3 and fig4 have their own entry scripts.
