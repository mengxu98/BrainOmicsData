#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

repo_dir <- normalizePath(value_after("--repo-dir", "."), mustWork = TRUE)
package_dir <- normalizePath(
  value_after("--package-dir", "../../data/BrainOmicsData/ScienceDB"),
  mustWork = TRUE
)
source_object <- normalizePath(
  value_after("--source-object", "../../data/BrainOmicsData/integration_25/objects_celltype_plot.rds"),
  mustWork = TRUE
)
source_metadata_arg <- value_after("--source-metadata", "")
assignment_file <- value_after("--celltype-assignment-file", "")
source_metadata_file <- if (nzchar(source_metadata_arg)) {
  normalizePath(source_metadata_arg, mustWork = TRUE)
} else {
  ""
}
source(file.path(repo_dir, "functions", "metadata_schema.R"))
source(file.path(repo_dir, "functions", "sample_schema.R"))
source(file.path(repo_dir, "sciencedb", "manifest.R"))
source(file.path(repo_dir, "sciencedb", "package_metadata.R"))
log_script <- file.path(repo_dir, "functions", "log_message.sh")


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

write_atomic_rds <- function(value, target) {
  tmp <- paste0(target, ".tmp.", Sys.getpid())
  raw_tmp <- paste0(tmp, ".uncompressed")
  on.exit(unlink(c(tmp, raw_tmp)), add = TRUE)
  pigz <- Sys.which("pigz")
  gzip <- Sys.which("gzip")
  if (nzchar(pigz) || nzchar(gzip)) {
    saveRDS(value, raw_tmp, compress = FALSE)
    compressor <- if (nzchar(pigz)) pigz else gzip
    threads <- suppressWarnings(as.integer(Sys.getenv(
      "BRAINOMICS_RDS_COMPRESSION_THREADS", "8"
    )))
    if (is.na(threads) || threads < 1L) threads <- 1L
    compressor_args <- if (nzchar(pigz)) {
      c("-n", "-p", as.character(threads), "-c", raw_tmp)
    } else {
      c("-n", "-c", raw_tmp)
    }
    status <- system2(compressor, compressor_args, stdout = tmp)
    if (!identical(status, 0L) || !file.exists(tmp) || file.size(tmp) == 0) {
      stop("Could not gzip the anonymized Seurat object")
    }
  } else {
    saveRDS(value, tmp, compress = TRUE)
  }
  if (!file.rename(tmp, target)) {
    stop("Could not atomically publish the anonymized Seurat object")
  }
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

canonicalize_missing_id <- function(values) {
  values <- as.character(values)
  missing <- is.na(values) | !nzchar(values) | values == "NA"
  values[missing] <- "NA"
  values
}

read_gzip_lines <- function(file) {
  con <- gzfile(file, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

stream_replace_first_column <- function(
  input, output, internal_ids, public_ids, chunk_size = 100000L
) {
  in_con <- gzfile(input, "rt")
  out_con <- gzfile(output, "wt")
  on.exit(close(in_con), add = TRUE)
  on.exit(close(out_con), add = TRUE)

  header <- readLines(in_con, n = 1L, warn = FALSE)
  if (length(header) != 1L) stop("empty file: ", input)
  writeLines(header, out_con)

  offset <- 0L
  mode <- NULL
  repeat {
    lines <- readLines(in_con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    tabs <- regexpr("\t", lines, fixed = TRUE)
    if (any(tabs < 1L)) stop("missing tab while processing: ", input)
    current <- substring(lines, 1L, tabs - 1L)
    rows <- (offset + 1L):(offset + length(lines))
    expected_internal <- internal_ids[rows]
    expected_public <- public_ids[rows]
    if (is.null(mode)) {
      mode <- if (identical(current, expected_internal)) {
        "internal"
      } else if (identical(current, expected_public)) {
        "public"
      } else {
        "invalid"
      }
    }
    expected <- if (identical(mode, "internal")) {
      expected_internal
    } else {
      expected_public
    }
    if (identical(mode, "invalid") || !identical(current, expected)) {
      bad <- which(current != expected)[1L]
      if (is.na(bad)) bad <- 1L
      stop(
        "cell order mismatch in ", input, " at data line ", offset + bad,
        ": observed=", current[[bad]], "; expected internal/public cell order"
      )
    }
    suffix <- substring(lines, tabs)
    writeLines(paste0(public_ids[rows], suffix), out_con)
    offset <- offset + length(lines)
  }
  if (offset != length(internal_ids)) {
    stop(
      "row count mismatch in ", input, ": observed=", offset,
      "; expected=", length(internal_ids)
    )
  }
}

validate_public_first_column <- function(
  input, public_ids, chunk_size = 100000L
) {
  in_con <- gzfile(input, "rt")
  on.exit(close(in_con), add = TRUE)
  header <- readLines(in_con, n = 1L, warn = FALSE)
  if (length(header) != 1L) stop("empty file: ", input)
  offset <- 0L
  repeat {
    lines <- readLines(in_con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    tabs <- regexpr("\t", lines, fixed = TRUE)
    if (any(tabs < 1L)) stop("missing tab while validating: ", input)
    current <- substring(lines, 1L, tabs - 1L)
    rows <- (offset + 1L):(offset + length(lines))
    if (!identical(current, public_ids[rows])) {
      stop("public cell order mismatch in ", input)
    }
    offset <- offset + length(lines)
  }
  if (offset != length(public_ids)) {
    stop(
      "row count mismatch in ", input, ": observed=", offset,
      "; expected=", length(public_ids)
    )
  }
  invisible(TRUE)
}

pair_map <- function(internal, public, label) {
  internal <- as.character(internal)
  public <- as.character(public)
  keep <- !is.na(internal) & internal != "NA" & nzchar(internal)
  pairs <- unique(data.frame(
    internal = internal[keep], public = public[keep],
    stringsAsFactors = FALSE
  ))
  if (anyNA(pairs$public) || any(pairs$public == "NA") ||
    any(!nzchar(pairs$public)) || anyDuplicated(pairs$internal) ||
    anyDuplicated(pairs$public)) {
    stop("Public ", label, " mapping is not one-to-one")
  }
  stats::setNames(pairs$public, pairs$internal)
}

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required")
})

metadata_file <- file.path(package_dir, "metadata", "metadata.tsv.gz")
shard_manifest_file <- file.path(
  package_dir, "expression", "shard_manifest.tsv"
)
integrated_pca_file <- file.path(package_dir, "embeddings", "integrated_pca.tsv.gz")
integrated_umap_file <- file.path(package_dir, "embeddings", "integrated_umap.tsv.gz")
unintegrated_umap_file <- file.path(package_dir, "embeddings", "unintegrated_umap.tsv.gz")
lisi_file <- file.path(package_dir, "validation", "lisi.tsv.gz")
object_file <- file.path(package_dir, "objects", "objects_celltype_plot.rds")
readme_file <- file.path(package_dir, "README.md")

thisutils::log_message("", "Loading public metadata")
metadata <- data.table::fread(metadata_file, sep = "\t", data.table = FALSE)
existing_metadata <- metadata
package_cells <- as.character(metadata$Cells)

thisutils::log_message("", "Loading source Seurat object")
object <- readRDS(source_object)
if (!inherits(object, "Seurat")) stop("source object is not a Seurat object")
object_cells <- colnames(object)
if (!identical(rownames(object@meta.data), object_cells)) {
  stop("Source object metadata and expression columns differ in cell order")
}
refreshed_source <- if (nzchar(source_metadata_file)) {
  thisutils::log_message("", "Loading refreshed source metadata")
  value <- readRDS(source_metadata_file)
  if (!is.data.frame(value) || anyDuplicated(rownames(value))) {
    stop("source metadata is not a uniquely row-named data frame")
  }
  source_rows <- rownames(value)
  value <- as.data.frame(value, stringsAsFactors = FALSE)
  rownames(value) <- source_rows
  value
} else {
  NULL
}
internal_cells <- if (is.null(refreshed_source)) {
  object_cells
} else {
  rownames(refreshed_source)
}
public_cells <- sprintf("Cell%07d", seq_along(internal_cells))
if (length(package_cells) != length(internal_cells) ||
  (!identical(object_cells, internal_cells) &&
    !identical(object_cells, public_cells))) {
  stop("Source object and package metadata differ in cell coverage or order")
}
package_is_internal <- identical(package_cells, internal_cells)
package_is_public <- identical(package_cells, public_cells)
if (!package_is_internal && !package_is_public) {
  stop("Package metadata contains neither internal nor canonical public cells")
}

source_metadata <- if (nzchar(source_metadata_file)) {
  refreshed_source[internal_cells, , drop = FALSE]
} else {
  object@meta.data[internal_cells, , drop = FALSE]
}
if (nzchar(source_metadata_file)) {
  refresh_columns <- c(
    "Cells", "Dataset", "Sample", "Original_Sample", "Original_Donor_ID",
    "Original_Sample_ID", "Original_Source_Sample_ID",
    "Donor_ID", "Global_Donor_ID", "Sample_ID", "Library_ID",
    "sequencing_technology", "Technology",
    "sequencing_modality_standardized", "Sequence", "Modality",
    "Age_num", "age_value", "Unit", "age_unit", "Age", "age_raw",
    "AgeIntervalID", "sex_standardized", "Sex", "BrainRegion",
    "Source_CellType", "source_cell_type_label", "CellType_raw",
    "Cluster", "Resolution_Cluster", "seurat_clusters", "CellType"
  )
  source_metadata <- source_metadata[
    , intersect(refresh_columns, names(source_metadata)),
    drop = FALSE
  ]
  source_metadata[] <- lapply(source_metadata, as.character)
  source_metadata <- add_metadata_schema(source_metadata)
  source_metadata <- add_sample_schema(source_metadata)
} else {
  source_metadata[] <- lapply(source_metadata, as.character)
  source_metadata <- add_metadata_schema(source_metadata)
  source_metadata <- add_sample_schema(source_metadata)
}
if (nzchar(assignment_file)) {
  annotation <- as.data.frame(readRDS(assignment_file))
  index <- match(internal_cells, annotation$Cells)
  stopifnot(nrow(annotation) == length(internal_cells), !anyNA(index), !anyDuplicated(annotation$Cells))
  source_metadata$Cluster <- as.character(annotation$Cluster[index])
  source_metadata$CellType <- as.character(annotation$CellType[index])
}
refreshed_metadata <- build_sciencedb_cell_metadata(
  source_metadata, internal_cells
)
preserved_fields <- c(
  "Dataset", "Original_Sample_ID", "AgeIntervalID", "Sex", "BrainRegion",
  "Source_CellType", "Cluster", "CellType"
)
if (nzchar(assignment_file)) preserved_fields <- setdiff(preserved_fields, c("Cluster", "CellType"))
identity_fields <- c(
  "Dataset", "Original_Sample_ID", "AgeIntervalID", "BrainRegion"
)
changed_identity_fields <- identity_fields[vapply(identity_fields, function(field) {
  !identical(
    canonicalize_missing_id(existing_metadata[[field]]),
    canonicalize_missing_id(refreshed_metadata[[field]])
  )
}, logical(1))]
if (package_is_public && length(changed_identity_fields)) {
  stop(
    "Refreshed source metadata changes fields outside the approved metadata refresh: ",
    paste(changed_identity_fields, collapse = ", ")
  )
}
if (package_is_public) {
  refreshed_metadata[preserved_fields] <- existing_metadata[preserved_fields]
}
required_id_columns <- c("Donor_ID", "Sample_ID", "Library_ID")
if (any(!required_id_columns %in% names(source_metadata)) ||
  any(!required_id_columns %in% names(existing_metadata))) {
  stop("Source/public metadata lack donor, sample or library identifiers")
}
if (package_is_internal) {
  donor_map <- make_id_map(source_metadata$Donor_ID, "D")
  sample_map <- make_id_map(source_metadata$Sample_ID, "S")
  library_map <- make_id_map(source_metadata$Library_ID, "L")
} else {
  donor_map <- pair_map(source_metadata$Donor_ID, existing_metadata$Donor_ID, "donor")
  sample_map <- pair_map(source_metadata$Sample_ID, existing_metadata$Sample_ID, "sample")
  library_map <- pair_map(source_metadata$Library_ID, existing_metadata$Library_ID, "library")
}
cell_map <- stats::setNames(public_cells, internal_cells)
metadata <- refreshed_metadata

shard_manifest <- data.table::fread(
  shard_manifest_file,
  sep = "\t", data.table = FALSE
)
required_shard_columns <- c(
  "Order", "Dataset", "Relative_Directory", "Cells", "Cell_Start",
  "Cell_End"
)
if (!all(required_shard_columns %in% names(shard_manifest)) ||
  !identical(as.integer(shard_manifest$Order), seq_len(nrow(shard_manifest))) ||
  sum(shard_manifest$Cells) != length(internal_cells) ||
  !identical(as.numeric(shard_manifest$Cell_Start), as.numeric(cumsum(c(
    1, head(shard_manifest$Cells, -1L)
  )))) ||
  !identical(
    as.numeric(shard_manifest$Cell_End),
    as.numeric(cumsum(shard_manifest$Cells))
  )) {
  stop("Expression shard manifest failed the public-cell order contract")
}
barcode_files <- file.path(
  package_dir, "expression", shard_manifest$Relative_Directory,
  "barcodes.tsv.gz"
)
if (any(!file.exists(barcode_files))) stop("Expression shard barcodes are missing")
barcodes <- unlist(lapply(barcode_files, read_gzip_lines), use.names = FALSE)
if (!identical(barcodes, internal_cells) &&
  !identical(barcodes, public_cells)) {
  stop("Expression shard barcodes contain neither internal nor public cell order")
}

thisutils::log_message("", "Anonymizing metadata and barcodes")
if (package_is_internal) {
  metadata$Cells <- public_cells
  metadata$Sample_ID <- replace_from_map(metadata$Sample_ID, sample_map)
  metadata$Donor_ID <- replace_from_map(metadata$Donor_ID, donor_map)
  metadata$Library_ID <- replace_from_map(metadata$Library_ID, library_map)
} else if (!identical(existing_metadata$Cells, public_cells) ||
  !identical(
    canonicalize_missing_id(existing_metadata$Donor_ID),
    canonicalize_missing_id(replace_from_map(
      source_metadata$Donor_ID, donor_map
    ))
  ) ||
  !identical(
    canonicalize_missing_id(existing_metadata$Sample_ID),
    canonicalize_missing_id(replace_from_map(
      source_metadata$Sample_ID, sample_map
    ))
  ) ||
  !identical(
    canonicalize_missing_id(existing_metadata$Library_ID),
    canonicalize_missing_id(replace_from_map(
      source_metadata$Library_ID, library_map
    ))
  )) {
  stop("Existing public metadata does not match the source-to-public ID map")
} else {
  metadata$Cells <- public_cells
  metadata$Sample_ID <- replace_from_map(source_metadata$Sample_ID, sample_map)
  metadata$Donor_ID <- replace_from_map(source_metadata$Donor_ID, donor_map)
  metadata$Library_ID <- replace_from_map(source_metadata$Library_ID, library_map)
}
write_atomic(metadata_file, function(tmp) write_tsv(metadata, tmp, gzip = TRUE))
for (i in seq_len(nrow(shard_manifest))) {
  rows <- seq.int(
    shard_manifest$Cell_Start[[i]], shard_manifest$Cell_End[[i]]
  )
  barcode_file <- barcode_files[[i]]
  write_atomic(barcode_file, function(tmp) {
    con <- gzfile(tmp, "wt")
    on.exit(close(con), add = TRUE)
    writeLines(public_cells[rows], con)
  })
}

for (file in c(integrated_pca_file, integrated_umap_file, unintegrated_umap_file, lisi_file)) {
  if (package_is_public) {
    thisutils::log_message("", paste("Validating existing public first-column cell IDs in", file))
    validate_public_first_column(file, public_cells)
  } else {
    thisutils::log_message("", paste("Anonymizing first-column cell IDs in", file))
    write_atomic(file, function(tmp) {
      stream_replace_first_column(file, tmp, internal_cells, public_cells)
    })
  }
}

thisutils::log_message("", "Renaming cells and metadata in Seurat object")
if (!identical(colnames(object), public_cells)) {
  object <- SeuratObject::RenameCells(object, new.names = public_cells)
}
metadata_for_object <- metadata[match(rownames(object@meta.data), metadata$Cells), , drop = FALSE]
if (anyNA(metadata_for_object$Cells)) stop("could not align anonymized metadata to Seurat object")
object@meta.data <- metadata_for_object
rownames(object@meta.data) <- metadata_for_object$Cells
if (ncol(object@meta.data) != 15L) {
  stop("The anonymized lightweight object must retain the 15-column public metadata schema")
}

thisutils::log_message("", "Saving anonymized Seurat object")
write_atomic_rds(object, object_file)

thisutils::log_message("", "Anonymizing reference identifiers in release sidecars")
table_files <- list.files(
  file.path(package_dir, c("metadata", "validation", "provenance")),
  recursive = TRUE, full.names = TRUE,
  pattern = "[.]tsv([.]gz)?$", ignore.case = TRUE
)
table_files <- setdiff(normalizePath(table_files, mustWork = TRUE), normalizePath(
  c(
    metadata_file, integrated_pca_file, integrated_umap_file,
    unintegrated_umap_file, lisi_file
  ),
  mustWork = TRUE
))
id_maps <- list(
  Global_Donor_ID = donor_map,
  Donor_ID = donor_map,
  Specimen_ID = sample_map,
  Sample_ID = sample_map,
  Library_ID = library_map
)
cell_columns <- c("Cells", "Cell", "cell_id")
table_dimensions <- vector("list", length(table_files))
for (table_file in table_files) {
  value <- data.table::fread(
    table_file,
    sep = "\t", data.table = FALSE,
    na.strings = c("", "NA")
  )
  table_dimensions[[match(table_file, table_files)]] <- data.frame(
    path = substring(table_file, nchar(package_dir) + 2L),
    rows = nrow(value),
    columns = ncol(value),
    nonzero_values = NA_real_,
    stringsAsFactors = FALSE
  )
  changed <- FALSE
  for (column in intersect(cell_columns, names(value))) {
    current <- as.character(value[[column]])
    nonmissing <- !is.na(current) & nzchar(current)
    internal_match <- nonmissing & current %in% names(cell_map)
    public_match <- nonmissing & current %in% unname(cell_map)
    if (any(internal_match)) {
      if (any(nonmissing & !internal_match & !public_match)) {
        stop("Mixed internal and unknown cell IDs in ", table_file)
      }
      current[internal_match] <- unname(cell_map[current[internal_match]])
      value[[column]] <- current
      changed <- TRUE
    }
  }
  for (column in intersect(names(id_maps), names(value))) {
    current <- as.character(value[[column]])
    mapping <- id_maps[[column]]
    matched <- !is.na(current) & current %in% names(mapping)
    if (any(matched)) {
      current[matched] <- unname(mapping[current[matched]])
      value[[column]] <- current
      changed <- TRUE
    }
  }
  file_columns <- grep("(^|_)File$", names(value), value = TRUE)
  for (column in file_columns) {
    current <- as.character(value[[column]])
    absolute <- !is.na(current) & startsWith(current, "/")
    if (any(absolute)) {
      current[absolute] <- basename(current[absolute])
      value[[column]] <- current
      changed <- TRUE
    }
  }
  if (changed) {
    write_atomic(table_file, function(tmp) {
      write_tsv(value, tmp, gzip = grepl("[.]gz$", table_file))
    })
  }
}
table_dimensions <- do.call(rbind, table_dimensions)

audit_index_file <- file.path(
  package_dir, "provenance", "audit_sidecar_index.tsv"
)
if (file.exists(audit_index_file)) unlink(audit_index_file)
audit_files <- list.files(
  file.path(package_dir, c("metadata", "validation", "provenance")),
  recursive = TRUE, full.names = TRUE
)
audit_files <- audit_files[!file.info(audit_files)$isdir]
# The manifest is regenerated after this index; it must not index itself.
audit_files <- audit_files[basename(audit_files) != "file_manifest.tsv"]
audit_files <- sort(normalizePath(audit_files, mustWork = TRUE))
audit_index <- data.frame(
  path = substring(audit_files, nchar(package_dir) + 2L),
  size_bytes = as.numeric(file.info(audit_files)$size),
  sha256 = unname(release_sha256(audit_files)),
  stringsAsFactors = FALSE
)
write_tsv(audit_index, audit_index_file)

thisutils::log_message("", "Updating README anonymization note")
readme <- readLines(readme_file, warn = FALSE)
note <- c(
  "",
  "Privacy and anonymization",
  "-------------------------",
  "The public package does not redistribute controlled-access raw human sequencing files. Project-internal cell, donor, biological-sample and library identifiers were replaced with package-internal anonymous identifiers, and that replacement map is not released. `Original_Sample_ID` retains only identifiers already published in the source repositories so that users can trace a released record back to its de-identified source record. Age and age-interval fields are retained because they are required scientific variables for age-stratified reuse."
)
if (!any(grepl("^Privacy and anonymization$", readme))) {
  writeLines(c(readme, note), readme_file)
}

thisutils::log_message("", "Updating file manifest and MD5 checksums")
old_manifest_file <- file.path(package_dir, "provenance", "file_manifest.tsv")
old_manifest <- if (file.exists(old_manifest_file)) {
  utils::read.delim(
    old_manifest_file,
    sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  data.frame()
}
dimension_audit <- if (all(
  c("path", "rows", "columns", "nonzero_values") %in% names(old_manifest)
)) {
  old_manifest[, c("path", "rows", "columns", "nonzero_values"), drop = FALSE]
} else {
  data.frame(
    path = character(), rows = numeric(), columns = numeric(),
    nonzero_values = numeric(), stringsAsFactors = FALSE
  )
}
dimension_audit <- rbind(dimension_audit, table_dimensions)
dimension_audit <- dimension_audit[
  !duplicated(dimension_audit$path, fromLast = TRUE), ,
  drop = FALSE
]
write_release_manifest(package_dir, dimension_audit)
thisutils::log_message("", "Public ScienceDB package anonymization completed")
