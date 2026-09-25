# Rasterize an assembled figure at a resolution that remains at least 600 dpi
# when Word places it at the specified final width.
export_pdf_png <- function(pdf, png, placed_width_in, min_dpi = 600L,
                           min_effective_dpi = min_dpi) {
  stopifnot(file.exists(pdf), placed_width_in > 0, min_dpi > 0,
            min_effective_dpi > 0)
  info <- system2("pdfinfo", pdf, stdout = TRUE)
  size <- grep("^Page size:", info, value = TRUE)
  if (length(size) != 1L) stop("Cannot read PDF page size: ", pdf)
  width_pt <- as.numeric(sub("^Page size:\\s*([0-9.]+).*", "\\1", size))
  dpi <- ceiling(max(min_dpi, min_effective_dpi * placed_width_in / (width_pt / 72)))
  prefix <- sub("\\.png$", "", png, ignore.case = TRUE)
  dir.create(dirname(png), recursive = TRUE, showWarnings = FALSE)
  status <- system2("pdftoppm", c("-f", "1", "-singlefile", "-png",
    "-r", as.character(dpi), shQuote(pdf), shQuote(prefix)))
  if (status != 0L || !file.exists(png)) stop("PNG export failed: ", png)
  # Export the same PDF at the same resolution, using lossless compression.
  tiff <- paste0(prefix, ".tiff")
  status <- system2("pdftoppm", c("-f", "1", "-singlefile", "-tiff",
    "-tiffcompression", "lzw", "-r", as.character(dpi),
    shQuote(pdf), shQuote(prefix)))
  if (status != 0L || !file.exists(paste0(prefix, ".tif")))
    stop("TIFF export failed: ", tiff)
  if (!file.rename(paste0(prefix, ".tif"), tiff)) stop("Cannot save TIFF: ", tiff)
  message(png, " exported at ", dpi, " dpi for ", placed_width_in, " in placement")
  message(tiff, " exported with LZW compression")
  invisible(png)
}
