#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

repo_dir <- normalizePath(value_after("--repo-dir", "."), mustWork = TRUE)
package_dir <- normalizePath(
  value_after("--package-dir", "../../data/BrainOmicsData/ScienceDB"),
  mustWork = TRUE
)
source_object <- normalizePath(
  value_after("--source-object", "../../data/BrainOmicsData/integration/objects_celltype_plot.rds"),
  mustWork = TRUE
)
log_script <- file.path(repo_dir, "functions", "log_message.sh")

log_message <- function(text, type = "info") {
  if (file.exists(log_script)) {
    status <- system2("bash", c(log_script, text, "--message-type", type), stdout = "", stderr = "")
    if (!identical(status, 0L)) message(text)
  } else {
    message(text)
  }
}

write_tsv <- function(x, file, gzip = grepl("\\.gz$", file)) {
  con <- if (gzip) gzfile(file, "wt") else file(file, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
}

write_atomic <- function(target, writer) {
  tmp <- paste0(target, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writer(tmp)
  file.rename(tmp, target)
}

make_id_map <- function(values, prefix) {
  values <- as.character(values)
  uniq <- unique(values[!is.na(values) & values != "NA" & nzchar(values)])
  width <- max(6, nchar(length(uniq)))
  ids <- sprintf("%s%0*d", prefix, width, seq_along(uniq))
  stats::setNames(ids, uniq)
}

replace_from_map <- function(values, map, missing_prefix = "NA") {
  values <- as.character(values)
  out <- unname(map[values])
  out[is.na(out) & (is.na(values) | values == "NA" | !nzchar(values))] <- missing_prefix
  out
}

read_gzip_lines <- function(file) {
  con <- gzfile(file, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

stream_replace_first_column <- function(input, output, old_ids, new_ids, chunk_size = 100000L) {
  in_con <- gzfile(input, "rt")
  out_con <- gzfile(output, "wt")
  on.exit(close(in_con), add = TRUE)
  on.exit(close(out_con), add = TRUE)

  header <- readLines(in_con, n = 1L, warn = FALSE)
  if (length(header) != 1L) stop("empty file: ", input)
  writeLines(header, out_con)

  offset <- 0L
  repeat {
    lines <- readLines(in_con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    tabs <- regexpr("\t", lines, fixed = TRUE)
    if (any(tabs < 1L)) stop("missing tab while processing: ", input)
    current_old <- substring(lines, 1L, tabs - 1L)
    expected <- old_ids[(offset + 1L):(offset + length(lines))]
    if (!identical(current_old, expected)) {
      bad <- which(current_old != expected)[1L]
      stop(
        "cell order mismatch in ", input, " at data line ", offset + bad,
        ": observed=", current_old[[bad]], "; expected=", expected[[bad]]
      )
    }
    suffix <- substring(lines, tabs)
    writeLines(paste0(new_ids[(offset + 1L):(offset + length(lines))], suffix), out_con)
    offset <- offset + length(lines)
  }
  if (offset != length(old_ids)) {
    stop("row count mismatch in ", input, ": observed=", offset, "; expected=", length(old_ids))
  }
}

update_manifest <- function(package_dir) {
  files <- list.files(package_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  files <- files[!basename(files) %in% c("file_manifest.tsv", "md5sum.txt")]
  rel <- sub(paste0("^", normalizePath(package_dir), "/?"), "", normalizePath(files))
  md5 <- unname(tools::md5sum(files))
  manifest <- data.frame(
    file = rel,
    size_bytes = file.info(files)$size,
    description = ifelse(grepl("^expression/", rel), "10X-compatible expression matrix component",
      ifelse(grepl("^metadata/", rel), "Harmonized anonymized metadata table",
        ifelse(grepl("^embeddings/", rel), "Anonymized embedding table",
          ifelse(grepl("^validation/", rel), "Anonymized integration validation table",
            ifelse(grepl("^objects/", rel), "Anonymized reusable Seurat object",
              ifelse(grepl("^provenance/", rel), "Source dataset references and citation provenance",
                ifelse(grepl("^scripts/", rel), "Example reader script", "Package documentation"))))))),
    md5 = md5,
    stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$file), ]
  write_tsv(manifest[, c("file", "size_bytes", "description")], file.path(package_dir, "file_manifest.tsv"), gzip = FALSE)
  writeLines(paste(manifest$md5, manifest$file), file.path(package_dir, "md5sum.txt"))
}

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required")
})

metadata_file <- file.path(package_dir, "metadata", "metadata.tsv.gz")
barcodes_file <- file.path(package_dir, "expression", "barcodes.tsv.gz")
integrated_pca_file <- file.path(package_dir, "embeddings", "integrated_pca.tsv.gz")
integrated_umap_file <- file.path(package_dir, "embeddings", "integrated_umap.tsv.gz")
unintegrated_umap_file <- file.path(package_dir, "embeddings", "unintegrated_umap.tsv.gz")
lisi_file <- file.path(package_dir, "validation", "lisi.tsv.gz")
object_file <- file.path(package_dir, "objects", "objects_celltype_plot.rds")
readme_file <- file.path(package_dir, "README.md")

log_message("Loading public metadata")
metadata <- data.table::fread(metadata_file, sep = "\t", data.table = FALSE)
old_cells <- metadata$Cells
old_barcodes <- read_gzip_lines(barcodes_file)
if (!identical(old_cells, old_barcodes)) {
  stop("metadata Cells and expression barcodes are not in the same order")
}

new_cells <- sprintf("Cell%07d", seq_along(old_cells))
cell_map <- stats::setNames(new_cells, old_cells)
donor_map <- make_id_map(metadata$Donor_ID, "D")
sample_map <- make_id_map(metadata$Sample_ID, "S")
library_map <- make_id_map(metadata$Library_ID, "L")

log_message("Anonymizing metadata and barcodes")
metadata$Cells <- new_cells
metadata$Sample <- replace_from_map(metadata$Sample, sample_map)
metadata$Sample_ID <- replace_from_map(metadata$Sample_ID, sample_map)
metadata$Donor_ID <- replace_from_map(metadata$Donor_ID, donor_map)
metadata$Library_ID <- replace_from_map(metadata$Library_ID, library_map)
write_atomic(metadata_file, function(tmp) write_tsv(metadata, tmp, gzip = TRUE))
write_atomic(barcodes_file, function(tmp) {
  con <- gzfile(tmp, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(new_cells, con)
})

for (file in c(integrated_pca_file, integrated_umap_file, unintegrated_umap_file, lisi_file)) {
  log_message(paste("Anonymizing first-column cell IDs in", file))
  write_atomic(file, function(tmp) stream_replace_first_column(file, tmp, old_cells, new_cells))
}

log_message("Loading source Seurat object")
object <- readRDS(source_object)
if (!inherits(object, "Seurat")) stop("source object is not a Seurat object")
if (!identical(colnames(object), old_cells)) {
  stop("source object cell order does not match public metadata")
}

log_message("Renaming cells and metadata in Seurat object")
object <- SeuratObject::RenameCells(object, new.names = new_cells)
metadata_for_object <- metadata[match(rownames(object@meta.data), metadata$Cells), , drop = FALSE]
if (anyNA(metadata_for_object$Cells)) stop("could not align anonymized metadata to Seurat object")
object@meta.data$Cells <- metadata_for_object$Cells
object@meta.data$Sample <- metadata_for_object$Sample
object@meta.data$Sample_ID <- metadata_for_object$Sample_ID
object@meta.data$Donor_ID <- metadata_for_object$Donor_ID
object@meta.data$Library_ID <- metadata_for_object$Library_ID
drop_cols <- intersect(c("Original_Sample", "Original_Sample_ID", "sample_schema_rule"), names(object@meta.data))
if (length(drop_cols)) object@meta.data <- object@meta.data[, setdiff(names(object@meta.data), drop_cols), drop = FALSE]

log_message("Saving anonymized Seurat object")
write_atomic(object_file, function(tmp) saveRDS(object, tmp, compress = TRUE))

log_message("Updating README anonymization note")
readme <- readLines(readme_file, warn = FALSE)
note <- c(
  "",
  "Privacy and anonymization",
  "-------------------------",
  "The public package does not redistribute controlled-access raw human sequencing files. Cell, donor, biological-sample and library identifiers in the released package were replaced with package-internal anonymous identifiers. The anonymization mapping is not released. Age and age-interval fields are retained because they are required scientific variables for age-stratified reuse."
)
if (!any(grepl("^Privacy and anonymization$", readme))) {
  writeLines(c(readme, note), readme_file)
}

log_message("Updating file manifest and MD5 checksums")
update_manifest(package_dir)
log_message("Public ScienceDB package anonymization completed")
