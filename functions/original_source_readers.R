require_original_source_packages <- function() {
  packages <- c("Matrix", "SeuratObject")
  missing <- packages[!vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]
  if (length(missing) > 0L) {
    stop("Missing original-source packages: ", paste(missing, collapse = ", "))
  }
}

required_source_files <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop("Required source files are missing: ", paste(missing, collapse = ", "))
  }
  invisible(files)
}

read_source_table <- function(file, sep = "\t") {
  required_source_files(file)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(
      data.table::fread(
        file,
        sep = sep,
        data.table = FALSE,
        check.names = FALSE,
        showProgress = FALSE
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  read.delim(
    file,
    sep = sep,
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_space_count_table_with_gene_column <- function(file) {
  required_source_files(file)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required for space-delimited count matrices")
  }
  connection <- if (grepl("[.]gz$", file, ignore.case = TRUE)) {
    gzfile(file, "rt")
  } else {
    file(file, "rt")
  }
  header <- tryCatch(
    readLines(connection, n = 1L, warn = FALSE),
    finally = close(connection)
  )
  cells <- strsplit(trimws(header), "[[:space:]]+")[[1L]]
  if (length(cells) == 0L || any(!nzchar(cells)) || anyDuplicated(cells)) {
    stop("Invalid cell header in space-delimited count matrix: ", file)
  }
  table <- data.table::fread(
    file,
    sep = " ",
    header = FALSE,
    skip = 1L,
    col.names = c("gene", cells),
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE
  )
  if (ncol(table) != length(cells) + 1L) {
    stop("Count matrix width differs from its cell header: ", file)
  }
  as.data.frame(
    table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_geo_family_soft_metadata <- function(file) {
  required_source_files(file)
  con <- if (grepl("[.]gz$", file, ignore.case = TRUE)) {
    gzfile(file, "rt")
  } else {
    file(file, "rt")
  }
  lines <- tryCatch(readLines(con, warn = FALSE), finally = close(con))
  starts <- which(startsWith(lines, "^SAMPLE = "))
  if (length(starts) == 0L) {
    stop("GEO family SOFT contains no sample blocks: ", file)
  }
  ends <- c(starts[-1L] - 1L, length(lines))
  rows <- Map(function(start, end) {
    block <- lines[start:end]
    accession <- sub("^\\^SAMPLE = ", "", block[[1L]])
    fields <- block[startsWith(block, "!Sample_")]
    split_field <- function(value) {
      pieces <- strsplit(value, " = ", fixed = TRUE)[[1L]]
      if (length(pieces) < 2L) {
        return(c(NA_character_, NA_character_))
      }
      c(sub("^!Sample_", "", pieces[[1L]]), paste(pieces[-1L], collapse = " = "))
    }
    parsed <- lapply(fields, split_field)
    keys <- vapply(parsed, `[[`, character(1), 1L)
    values <- vapply(parsed, `[[`, character(1), 2L)
    ordinary <- keys != "characteristics_ch1"
    row <- list(Sample_geo_accession = accession)
    for (index in which(ordinary)) {
      key <- keys[[index]]
      target <- if (key %in% names(row)) {
        paste0(key, "_", sum(names(row) == key) + 1L)
      } else {
        key
      }
      row[[target]] <- values[[index]]
    }
    for (value in values[!ordinary]) {
      key <- trimws(sub(":.*$", "", value))
      target <- paste0(
        "characteristic_",
        gsub("[^A-Za-z0-9]+", "_", tolower(key))
      )
      if (target %in% names(row)) {
        stop("Duplicated GEO characteristic in ", accession, ": ", target)
      }
      row[[target]] <- trimws(sub("^[^:]+:[[:space:]]*", "", value))
    }
    as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  }, starts, ends)
  bind_rows_by_name(rows)
}

read_one_column_csv <- function(file) {
  table <- read.csv(
    gzfile(file),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(table) < 2L) {
    stop("Expected an index and value column in: ", file)
  }
  as.character(table[[ncol(table)]])
}

read_10x_triplet <- function(matrix_file, barcode_file, feature_file) {
  required_source_files(c(matrix_file, barcode_file, feature_file))
  counts <- Matrix::readMM(gzfile(matrix_file))
  counts <- methods::as(counts, "CsparseMatrix")
  if (!inherits(counts, "dgCMatrix")) {
    stop("10x Matrix Market input did not produce a dgCMatrix: ", matrix_file)
  }
  barcodes <- readLines(gzfile(barcode_file), warn = FALSE)
  features <- read.delim(
    gzfile(feature_file),
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  feature_names <- as.character(features[[min(2L, ncol(features))]])
  feature_names <- make.unique(feature_names)
  if (nrow(counts) != length(feature_names) ||
    ncol(counts) != length(barcodes)) {
    stop("10x matrix dimensions differ from feature/barcode files")
  }
  rownames(counts) <- feature_names
  colnames(counts) <- barcodes
  list(counts = counts, barcodes = barcodes, features = features)
}

discover_10x_triplets <- function(raw_dir) {
  matrices <- sort(list.files(
    raw_dir,
    pattern = "_matrix[.]mtx[.]gz$",
    full.names = TRUE
  ))
  if (length(matrices) == 0L) {
    stop("No extracted 10x matrices were found under ", raw_dir)
  }
  lapply(matrices, function(matrix_file) {
    stem <- sub("_matrix[.]mtx[.]gz$", "", matrix_file)
    barcodes <- paste0(stem, "_barcodes.tsv.gz")
    feature_candidates <- c(
      paste0(stem, "_features.tsv.gz"),
      paste0(stem, "_genes.tsv.gz")
    )
    feature_file <- feature_candidates[file.exists(feature_candidates)][1L]
    if (is.na(feature_file)) {
      stop("No feature/gene file found for ", matrix_file)
    }
    list(
      stem = basename(stem),
      matrix = matrix_file,
      barcodes = barcodes,
      features = feature_file
    )
  })
}

align_layer_features <- function(layers) {
  if (length(layers) == 0L || any(!vapply(
    layers,
    inherits,
    logical(1),
    what = "dgCMatrix"
  ))) {
    stop("Count layers must be a non-empty list of dgCMatrix objects")
  }
  features <- unique(unlist(lapply(layers, rownames), use.names = FALSE))
  aligned <- lapply(layers, function(layer) {
    if (identical(rownames(layer), features)) {
      return(layer)
    }
    missing <- setdiff(features, rownames(layer))
    if (length(missing) > 0L) {
      zero <- Matrix::sparseMatrix(
        i = integer(),
        j = integer(),
        dims = c(length(missing), ncol(layer)),
        dimnames = list(missing, colnames(layer))
      )
      layer <- rbind(layer, zero)
    }
    layer[features, , drop = FALSE]
  })
  names(aligned) <- names(layers)
  aligned
}

create_complete_source_object <- function(layers, metadata, dataset) {
  require_original_source_packages()
  if (is.null(names(layers)) || any(names(layers) == "") ||
    anyDuplicated(names(layers))) {
    stop("Source layers must have unique non-empty names")
  }
  layers <- align_layer_features(layers)
  cells <- unlist(lapply(layers, colnames), use.names = FALSE)
  if (anyDuplicated(cells) || anyNA(cells)) {
    stop("Source count layers contain duplicated or missing cell IDs")
  }
  if (!"Cells" %in% names(metadata) || anyDuplicated(metadata$Cells) ||
    !setequal(metadata$Cells, cells)) {
    stop("Source metadata does not cover every matrix cell exactly once")
  }
  metadata <- metadata[match(cells, metadata$Cells), , drop = FALSE]
  rownames(metadata) <- metadata$Cells
  object <- if (length(layers) == 1L) {
    SeuratObject::CreateSeuratObject(
      counts = layers[[1L]],
      meta.data = metadata,
      project = dataset
    )
  } else {
    SeuratObject::CreateSeuratObject(
      counts = layers,
      meta.data = metadata,
      project = dataset
    )
  }
  if (!setequal(colnames(object), cells) || ncol(object) != length(cells)) {
    stop("Seurat object changed complete source cell membership")
  }
  object
}

read_hca_submission_sheet <- function(file, sheet) {
  required_source_files(file)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("readxl is required to read original HCA metadata: ", file)
  }
  raw <- as.data.frame(
    readxl::read_excel(
      file,
      sheet = sheet,
      col_names = FALSE,
      .name_repair = "minimal"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  schema_row <- which(apply(raw, 1L, function(row) {
    any(startsWith(as.character(row), paste0(tolower(gsub(" ", "_", sheet)), ".")))
  }))
  if (length(schema_row) != 1L) {
    stop("Cannot identify the HCA schema row in sheet: ", sheet)
  }
  marker_row <- which(apply(raw, 1L, function(row) {
    any(as.character(row) == "FILL OUT INFORMATION BELOW THIS ROW", na.rm = TRUE)
  }))
  if (length(marker_row) != 1L || marker_row <= schema_row) {
    stop("Cannot identify the HCA data marker in sheet: ", sheet)
  }
  column_names <- as.character(unlist(raw[schema_row, ], use.names = FALSE))
  missing_name <- is.na(column_names) | !nzchar(column_names)
  column_names[missing_name] <- paste0("unused_", which(missing_name))
  column_names <- make.unique(column_names)
  result <- raw[seq.int(marker_row + 1L, nrow(raw)), , drop = FALSE]
  names(result) <- column_names
  keep <- apply(result, 1L, function(row) {
    any(!is.na(row) & nzchar(trimws(as.character(row))))
  })
  result[keep, , drop = FALSE]
}

normalize_gse202210_donor_key <- function(value) {
  key <- toupper(trimws(as.character(value)))
  key <- sub("_HHT$", "", key)
  gsub("[^A-Z0-9]", "", key)
}

read_gse202210_hca_donors <- function(file) {
  donors <- read_hca_submission_sheet(file, "Donor organism")
  required <- c(
    "donor_organism.biomaterial_core.biomaterial_id",
    "donor_organism.sex",
    "donor_organism.organism_age",
    "donor_organism.organism_age_unit.text",
    "donor_organism.human_specific.ethnicity.text",
    "donor_organism.diseases.text"
  )
  missing <- setdiff(required, names(donors))
  if (length(missing) > 0L) {
    stop("GSE202210 HCA donor metadata is missing: ", paste(missing, collapse = ", "))
  }
  donors$HCA_Donor_ID <- as.character(
    donors[["donor_organism.biomaterial_core.biomaterial_id"]]
  )
  donors$Donor_Key <- normalize_gse202210_donor_key(donors$HCA_Donor_ID)
  donors$Age_Years <- suppressWarnings(as.numeric(
    donors[["donor_organism.organism_age"]]
  ))
  donors$Age_Unit <- tolower(as.character(
    donors[["donor_organism.organism_age_unit.text"]]
  ))
  donors$Sex_Original <- as.character(donors[["donor_organism.sex"]])
  donors$Ethnicity_Original <- as.character(
    donors[["donor_organism.human_specific.ethnicity.text"]]
  )
  donors$Diagnosis_Original <- as.character(
    donors[["donor_organism.diseases.text"]]
  )
  if (anyNA(donors$Donor_Key) || any(!nzchar(donors$Donor_Key)) ||
    anyDuplicated(donors$Donor_Key)) {
    stop("GSE202210 HCA donor identifiers are incomplete or duplicated")
  }
  if (any(is.na(donors$Age_Years)) || any(donors$Age_Unit != "year")) {
    stop("GSE202210 HCA donor ages are incomplete or not reported in years")
  }
  donors[, c(
    "HCA_Donor_ID", "Donor_Key", "Age_Years", "Age_Unit",
    "Sex_Original", "Ethnicity_Original", "Diagnosis_Original"
  ), drop = FALSE]
}

write_source_reconstruction_audit <- function(object, output_file, dataset) {
  layers <- processed_count_layers(object)
  audit <- data.frame(
    Dataset = dataset,
    Complete_Cells = ncol(object),
    Complete_Features = nrow(object),
    Count_Layers = length(layers),
    Matrix_Class = paste(unique(vapply(
      layers,
      function(x) class(x)[[1L]],
      character(1)
    )), collapse = ";"),
    Cell_IDs_Unique = !anyDuplicated(colnames(object)),
    Metadata_Aligned = identical(colnames(object), rownames(object[[]])),
    All_Counts_Nonnegative = all(vapply(
      layers,
      function(x) all(x@x >= 0 & is.finite(x@x)),
      logical(1)
    )),
    stringsAsFactors = FALSE
  )
  write.table(
    audit,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  invisible(audit)
}
