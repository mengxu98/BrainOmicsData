# Rasterize an assembled figure at a resolution that remains at least 600 dpi
# when Word places it at the specified final width.
export_pdf_png <- function(pdf, png, placed_width_in, min_dpi = 600L) {
  stopifnot(file.exists(pdf), placed_width_in > 0, min_dpi > 0)
  info <- system2("pdfinfo", pdf, stdout = TRUE)
  size <- grep("^Page size:", info, value = TRUE)
  if (length(size) != 1L) stop("Cannot read PDF page size: ", pdf)
  width_pt <- as.numeric(sub("^Page size:\\s*([0-9.]+).*", "\\1", size))
  dpi <- ceiling(max(min_dpi, min_dpi * placed_width_in / (width_pt / 72)))
  prefix <- sub("\\.png$", "", png, ignore.case = TRUE)
  dir.create(dirname(png), recursive = TRUE, showWarnings = FALSE)
  status <- system2("pdftoppm", c("-f", "1", "-singlefile", "-png",
    "-r", as.character(dpi), shQuote(pdf), shQuote(prefix)))
  if (status != 0L || !file.exists(png)) stop("PNG export failed: ", png)
  message(png, " exported at ", dpi, " dpi for ", placed_width_in, " in placement")
  invisible(png)
}
