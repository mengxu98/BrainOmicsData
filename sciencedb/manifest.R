release_manifest_schema_version <- function() "2.0.0"

release_file_description <- function(path) {
  ifelse(grepl("^expression/shard_manifest[.]tsv$", path),
    "Ordered complete-expression shard inventory",
    ifelse(grepl("^expression/", path), "Dataset-resolved 10X-compatible expression shard component",
      ifelse(grepl("^metadata/", path), "Harmonized metadata or metadata audit table",
        ifelse(grepl("^embeddings/", path), "Cell-ordered embedding table",
          ifelse(grepl("^validation/", path), "Integration or annotation validation result",
            ifelse(grepl("^objects/", path), "Reusable lightweight Seurat object",
              ifelse(grepl("^provenance/", path), "Source reference, access and provenance record",
                ifelse(grepl("^scripts/", path), "Example reader script", "Package documentation")
              )
            )
          )
        )
      )
    )
  )
}

release_file_license_scope <- function(path) {
  ifelse(
    grepl("^scripts/", path),
    "MIT License",
    "Derived data: source-specific access and reuse terms in provenance/dataset_manifest.tsv"
  )
}

release_file_schema <- function(path) {
  ifelse(
    path == "metadata/metadata.tsv.gz",
    paste0("BrainOmicsData metadata schema ", brainomics_metadata_schema_version()),
    ifelse(
      grepl("^metadata/", path),
      "BrainOmicsData metadata sidecar schema 1.0.0",
      ifelse(
        grepl("^embeddings/", path),
        "BrainOmicsData cell-embedding schema 1.0.0",
        ifelse(
          grepl("^validation/", path),
          "BrainOmicsData validation schema 1.0.0",
          ifelse(
            grepl("^expression/", path),
            "BrainOmicsData dataset-resolved Matrix Market shard schema 2.0.0",
            "not applicable"
          )
        )
      )
    )
  )
}

release_sha256 <- function(files) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest is required to generate the release SHA-256 chain")
  }
  vapply(
    files,
    digest::digest,
    character(1L),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
}

write_release_manifest <- function(package_dir, dimensions = NULL) {
  package_dir <- normalizePath(package_dir, mustWork = TRUE)
  excluded <- c("file_manifest.tsv", "provenance/file_manifest.tsv", "md5sum.txt", "sha256sum.txt")
  files <- list.files(
    package_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE,
    no.. = TRUE
  )
  files <- files[!file.info(files)$isdir]
  files <- files[!basename(files) %in% excluded]
  files <- sort(normalizePath(files, mustWork = TRUE))
  prefix_length <- nchar(package_dir) + 2L
  paths <- substring(files, prefix_length)

  manifest <- data.frame(
    path = paths,
    file_name = basename(paths),
    size_bytes = as.numeric(file.info(files)$size),
    rows = NA_real_,
    columns = NA_real_,
    nonzero_values = NA_real_,
    schema_version = release_file_schema(paths),
    license_scope = release_file_license_scope(paths),
    description = release_file_description(paths),
    md5 = unname(tools::md5sum(files)),
    sha256 = unname(release_sha256(files)),
    stringsAsFactors = FALSE
  )
  if (!is.null(dimensions) && nrow(dimensions) > 0L) {
    required <- c("path", "rows", "columns", "nonzero_values")
    if (!all(required %in% names(dimensions)) || anyDuplicated(dimensions$path)) {
      stop("Release dimension audit must contain one row per path")
    }
    matched <- match(manifest$path, dimensions$path)
    for (column in required[-1L]) {
      use <- !is.na(matched)
      manifest[[column]][use] <- dimensions[[column]][matched[use]]
    }
  }

  manifest_file <- file.path(package_dir, "provenance", "file_manifest.tsv")
  dir.create(file.path(package_dir, "provenance"), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    manifest,
    manifest_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
    na = "NA"
  )
  checksum_paths <- c(paths, "provenance/file_manifest.tsv")
  manifest_md5 <- unname(tools::md5sum(manifest_file))
  manifest_sha256 <- unname(release_sha256(manifest_file))
  writeLines(
    paste(c(manifest$md5, manifest_md5), checksum_paths),
    file.path(package_dir, "md5sum.txt")
  )
  invisible(manifest)
}

write_release_manifest_incremental <- function(
  package_dir, previous_manifest, dimensions = NULL,
  changed_paths = character()
) {
  package_dir <- normalizePath(package_dir, mustWork = TRUE)
  required_previous <- c(
    "path", "file_name", "size_bytes", "rows", "columns",
    "nonzero_values", "schema_version", "license_scope", "description",
    "md5", "sha256"
  )
  if (!identical(names(previous_manifest), required_previous) ||
    anyDuplicated(previous_manifest$path) ||
    anyNA(previous_manifest[, c("path", "md5", "sha256")])) {
    stop("Previous release manifest cannot support an incremental refresh")
  }
  excluded <- c("file_manifest.tsv", "provenance/file_manifest.tsv", "md5sum.txt", "sha256sum.txt")
  files <- list.files(
    package_dir,
    recursive = TRUE, full.names = TRUE,
    all.files = FALSE, no.. = TRUE
  )
  files <- files[!file.info(files)$isdir]
  files <- files[!basename(files) %in% excluded]
  files <- sort(normalizePath(files, mustWork = TRUE))
  paths <- substring(files, nchar(package_dir) + 2L)
  sizes <- as.numeric(file.info(files)$size)
  previous_index <- match(paths, previous_manifest$path)
  reuse <- !is.na(previous_index) &
    as.numeric(previous_manifest$size_bytes[previous_index]) == sizes &
    !paths %in% changed_paths

  md5 <- sha256 <- rep(NA_character_, length(files))
  md5[reuse] <- previous_manifest$md5[previous_index[reuse]]
  sha256[reuse] <- previous_manifest$sha256[previous_index[reuse]]
  refresh <- !reuse
  if (any(refresh)) {
    md5[refresh] <- unname(tools::md5sum(files[refresh]))
    sha256[refresh] <- unname(release_sha256(files[refresh]))
  }
  manifest <- data.frame(
    path = paths,
    file_name = basename(paths),
    size_bytes = sizes,
    rows = NA_real_, columns = NA_real_, nonzero_values = NA_real_,
    schema_version = release_file_schema(paths),
    license_scope = release_file_license_scope(paths),
    description = release_file_description(paths),
    md5 = md5, sha256 = sha256,
    stringsAsFactors = FALSE
  )
  if (!is.null(dimensions) && nrow(dimensions) > 0L) {
    required <- c("path", "rows", "columns", "nonzero_values")
    if (!all(required %in% names(dimensions)) ||
      anyDuplicated(dimensions$path)) {
      stop("Release dimension audit must contain one row per path")
    }
    matched <- match(manifest$path, dimensions$path)
    for (column in required[-1L]) {
      use <- !is.na(matched)
      manifest[[column]][use] <- dimensions[[column]][matched[use]]
    }
  }

  manifest_file <- file.path(package_dir, "provenance", "file_manifest.tsv")
  dir.create(file.path(package_dir, "provenance"), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    manifest, manifest_file,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  checksum_paths <- c(paths, "provenance/file_manifest.tsv")
  manifest_md5 <- unname(tools::md5sum(manifest_file))
  manifest_sha256 <- unname(release_sha256(manifest_file))
  writeLines(
    paste(c(manifest$md5, manifest_md5), checksum_paths),
    file.path(package_dir, "md5sum.txt")
  )
  message(
    "Incremental manifest refresh reused ", sum(reuse),
    " payload hashes and recomputed ", sum(refresh)
  )
  invisible(manifest)
}
