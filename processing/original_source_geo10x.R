process_geo_10x <- function(dataset_name) {
  triplets <- discover_10x_triplets(raw_dir)
  soft_file <- file.path(raw_dir, paste0(dataset_name, "_family.soft.gz"))
  geo <- read_geo_family_soft_metadata(soft_file)
  layers <- list()
  rows <- list()
  for (record in triplets) {
    if (dataset_name == "GSE178175" &&
      !grepl("_10X_sc$", record$stem)) {
      next
    }
    triplet <- read_10x_triplet(
      record$matrix,
      record$barcodes,
      record$features
    )
    gsm <- sub("_.*$", "", record$stem)
    geo_index <- match(gsm, geo$Sample_geo_accession)
    if (is.na(geo_index)) stop("No GEO metadata for ", gsm)
    geo_row <- geo[geo_index, , drop = FALSE]
    sample <- if (dataset_name == "GSE178175") {
      sub("_10X_sc$", "", sub("^[^_]+_", "", record$stem))
    } else {
      sub("^[^_]+_", "", record$stem)
    }
    cells <- paste(sample, triplet$barcodes, sep = "_")
    colnames(triplet$counts) <- cells
    diagnosis_column <- if (dataset_name == "GSE178175") {
      "characteristic_disease_state"
    } else {
      "characteristic_diagnosis"
    }
    diagnosis <- rep(as.character(geo_row[[diagnosis_column]]), length(cells))
    region <- if (dataset_name == "GSE178175") {
      "Frontal cortex"
    } else {
      "Prefrontal cortex"
    }
    age <- if (dataset_name == "GSE178175") {
      paste0(geo_row$characteristic_age, " years")
    } else {
      rep(NA_character_, length(cells))
    }
    sex <- if (dataset_name == "GSE178175") {
      rep(geo_row$characteristic_sex, length(cells))
    } else {
      rep(NA_character_, length(cells))
    }
    meta <- base_source_metadata(
      cells,
      sample,
      sample,
      gsm,
      gsm,
      age,
      sex,
      region,
      diagnosis,
      NA_character_
    )
    meta$Source_Barcode <- triplet$barcodes
    layers[[sample]] <- triplet$counts
    rows[[sample]] <- meta
  }
  if (length(layers) == 0L) stop("No eligible 10x RNA matrices were found")
  layers <- align_layer_features(layers)
  metadata <- bind_rows_by_name(rows)
  metadata$Source_CellType_Original <- NA_character_
  metadata$CellType_raw <- NA_character_
  metadata$Source_CellType_Annotation_Origin <-
    "not_released_by_original_study"
  metadata$Library_Chemistry <- "10x Chromium 3' v3"
  metadata$Sequencing_Platform <- "Illumina NovaSeq 6000"
  if (dataset_name == "GSE178175") {
    metadata$Analysis_Include <- TRUE
    metadata$Analysis_Role <- "reference"
    metadata$Exclusion_Reason <- NA_character_
    metadata <- age_provenance(
      metadata,
      metadata$Age,
      sub(" years$", "", metadata$Age),
      "years",
      "postnatal age at death",
      "GEO GSE178175 family SOFT",
      "GEO age retained in years",
      "source GEO sample metadata",
      FALSE
    )
  } else {
    hca_file <- file.path(
      raw_dir,
      "ZhangLabPdBrainNuclei_metadata_10-01-2024.xlsx"
    )
    hca_donors <- read_gse202210_hca_donors(hca_file)
    geo_donor_id <- metadata$Original_Donor_ID
    donor_index <- match(
      normalize_gse202210_donor_key(geo_donor_id),
      hca_donors$Donor_Key
    )
    if (anyNA(donor_index)) {
      stop("A GSE202210 GEO donor does not map to the original HCA metadata")
    }
    if (length(unique(donor_index)) != nrow(hca_donors) ||
      nrow(hca_donors) != 12L) {
      stop("GSE202210 GEO/HCA donor mapping is not one-to-one for 12 donors")
    }
    geo_control <- grepl("control", metadata$Diagnosis, ignore.case = TRUE)
    hca_control <- grepl(
      "normal",
      hca_donors$Diagnosis_Original[donor_index],
      ignore.case = TRUE
    )
    if (!identical(geo_control, hca_control)) {
      stop("GSE202210 GEO and original HCA diagnosis fields disagree")
    }
    metadata$Original_Donor_ID_GEO <- geo_donor_id
    metadata$Original_Donor_ID <- hca_donors$HCA_Donor_ID[donor_index]
    metadata$Ethnicity_raw <- hca_donors$Ethnicity_Original[donor_index]
    metadata$Age <- paste0(
      format_canonical_age_number(hca_donors$Age_Years[donor_index]),
      " years"
    )
    source_sex <- tolower(hca_donors$Sex_Original[donor_index])
    metadata$Sex <- ifelse(
      source_sex == "female",
      "Female",
      ifelse(source_sex == "male", "Male", NA_character_)
    )
    if (anyNA(metadata$Sex)) {
      stop("GSE202210 original HCA donor sex is incomplete")
    }
    metadata$Analysis_Include <- grepl(
      "control",
      metadata$Diagnosis,
      ignore.case = TRUE
    )
    metadata$Analysis_Role <- ifelse(
      metadata$Analysis_Include,
      "reference",
      "reference_excluded"
    )
    metadata$Exclusion_Reason <- ifelse(
      metadata$Analysis_Include,
      NA_character_,
      "Parkinson disease source donor"
    )
    metadata <- age_provenance(
      metadata,
      metadata$Age,
      as.character(hca_donors$Age_Years[donor_index]),
      "years",
      "postnatal age at death",
      paste(
        "original HCA project submission",
        "9a23ac2d-93dd-4bac-9bb8-040e6426db9d;",
        "Zhu et al. 2024, DOI 10.1126/scitranslmed.abo1997"
      ),
      "original HCA age retained in years",
      "high: original project donor metadata",
      FALSE
    )
  }
  create_complete_source_object(layers, metadata, dataset_name)
}
