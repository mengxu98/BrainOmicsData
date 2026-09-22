source("functions/utils.R")

existing_reference_datasets <- function() {
  c(
    "AllenM1", "EGAS00001006537",
    "GSE104276", "GSE168408",
    "GSE186538",
    "GSE204683", "GSE207334", "GSE212606", "GSE217511",
    "GSE296073", "GSE67835", "GSE81475",
    "GSE97942", "HYPOMAP", "Li_et_al_2018", "Ma_et_al_2022",
    "PRJCA015229", "ROSMAP", "SomaMut"
  )
}

excluded_reference_datasets <- function() {
  c(
    "GSE103723", "GSE199762", "GSE261983",
    "Nowakowski_et_al_2017"
  )
}

additional_reference_datasets <- function() {
  c("GSE294786", "Velmeshev_2023", "Wang_2025")
}

reference_datasets <- function() {
  c(existing_reference_datasets(), additional_reference_datasets())
}

query_validation_datasets <- function() {
  "EGAD00001006049"
}

formal_datasets <- function() {
  c(reference_datasets(), query_validation_datasets())
}

configure_complete_integration_future <- function() {
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("The future package is required for complete-reference integration")
  }
  future::plan(future::sequential)
  workers <- future::nbrOfWorkers()
  if (!identical(as.integer(workers), 1L)) {
    stop(
      "Complete-reference integration must use one in-process future worker"
    )
  }
  # Seurat's RPCA integration uses future.apply even with a sequential plan.
  # The complete reference legitimately closes over objects larger than the
  # default globals limit, but sequential execution does not transfer them to
  # another R process. Disable that transfer guard only after confirming that
  # execution is in-process and single-worker.
  options(future.globals.maxSize = Inf)
  list(
    Plan = "sequential",
    Workers = workers,
    Globals_Max_Size = "Inf (in-process; no worker transfer)"
  )
}

assign_integration_batch <- function(
  meta,
  model = c("dataset", "verified_technical")
) {
  model <- match.arg(model)
  required <- c(
    "Dataset", "Global_Donor_ID", "Specimen_ID",
    "Technical_Batch_ID", "Technical_Batch_Verification_Status",
    "Technical_Batch_Evidence", "Technical_Batch_Source_Column",
    "Technical_Batch_Derivation_Rule", "Technical_Batch_Use_Status",
    "Technical_Batch_Use_Evidence"
  )
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L) {
    stop(
      "integration metadata is missing audited batch columns: ",
      paste(missing, collapse = ", ")
    )
  }
  dataset <- normalize_missing_metadata(meta$Dataset)
  if (anyNA(dataset)) {
    stop("Dataset is missing for one or more integration cells")
  }
  dataset_key <- paste0("dataset:", dataset)
  technical_batch <- normalize_missing_metadata(meta$Technical_Batch_ID)
  technical_status <- tolower(as.character(
    meta$Technical_Batch_Verification_Status
  ))
  technical_evidence <- normalize_missing_metadata(
    meta$Technical_Batch_Evidence
  )
  technical_use_status <- tolower(as.character(
    meta$Technical_Batch_Use_Status
  ))
  technical_use_evidence <- normalize_missing_metadata(
    meta$Technical_Batch_Use_Evidence
  )
  verified <- !is.na(technical_batch) &
    grepl("^verified([ :]|$)", technical_status) &
    !is.na(technical_evidence)
  curated_eligible <- verified &
    grepl("^eligible([ :]|$)", technical_use_status) &
    !is.na(technical_use_evidence)

  biological_unit <- normalize_missing_metadata(meta$Global_Donor_ID)
  missing_donor <- is.na(biological_unit)
  biological_unit[missing_donor] <- normalize_missing_metadata(
    meta$Specimen_ID[missing_donor]
  )
  batch_units <- stats::setNames(integer(), character())
  if (any(verified)) {
    split_units <- split(
      biological_unit[verified],
      technical_batch[verified]
    )
    batch_units <- vapply(
      split_units,
      function(value) length(unique(stats::na.omit(value))),
      integer(1)
    )
  }
  spans_multiple_biological_units <- verified &
    unname(batch_units[technical_batch]) >= 2L
  spans_multiple_biological_units[is.na(
    spans_multiple_biological_units
  )] <- FALSE

  dataset_eligible <- rep(FALSE, nrow(meta))
  dataset_eligibility_reason <- rep(
    "dataset fallback: no paper/repository-verified technical factor",
    nrow(meta)
  )
  for (dataset_value in unique(dataset)) {
    index <- which(dataset == dataset_value)
    complete_verified <- all(verified[index])
    complete_curated_eligibility <- all(curated_eligible[index])
    technical_levels <- unique(stats::na.omit(
      technical_batch[index]
    ))
    multiple_levels <- length(technical_levels) >= 2L
    all_levels_span_units <- complete_verified &&
      complete_curated_eligibility && multiple_levels &&
      all(spans_multiple_biological_units[index])
    dataset_eligible[index] <- all_levels_span_units
    dataset_eligibility_reason[index] <- if (all_levels_span_units) {
      paste(
        "eligible: complete paper/repository-verified technical factor;",
        "at least two levels; every level spans at least two biological units"
      )
    } else if (!any(verified[index])) {
      "dataset fallback: no paper/repository-verified technical factor"
    } else if (!complete_verified) {
      "dataset fallback: verified technical factor is incomplete"
    } else if (!complete_curated_eligibility) {
      paste(
        "dataset fallback: verified factor is recorded but not eligible",
        "for correction because of the prespecified biological-confounding",
        "or design audit"
      )
    } else if (!multiple_levels) {
      "dataset fallback: fewer than two technical-batch levels"
    } else {
      paste(
        "dataset fallback: at least one technical level is confounded",
        "with a single donor/specimen"
      )
    }
  }
  use_technical <- model == "verified_technical" & dataset_eligible
  technical_model_key <- dataset_key
  technical_model_key[dataset_eligible] <- technical_batch[dataset_eligible]
  meta$Integration_Batch_ID_Dataset <- dataset_key
  meta$Integration_Batch_ID_Verified_Technical <- technical_model_key
  meta$Integration_Batch_ID <- if (model == "dataset") {
    meta$Integration_Batch_ID_Dataset
  } else {
    meta$Integration_Batch_ID_Verified_Technical
  }
  meta$Integration_Batch_Model <- model
  meta$Integration_Batch_Source <- rep("source dataset", nrow(meta))
  if (model == "verified_technical") {
    meta$Integration_Batch_Source <- dataset_eligibility_reason
    meta$Integration_Batch_Source[use_technical] <-
      "paper/repository-verified technical batch"
  }
  meta$Technical_Batch_Spans_Multiple_Biological_Units <-
    spans_multiple_biological_units
  meta$Technical_Batch_Dataset_Eligible <- dataset_eligible
  meta$Technical_Batch_Dataset_Eligibility_Reason <-
    dataset_eligibility_reason
  if (anyNA(meta$Integration_Batch_ID) ||
    any(meta$Integration_Batch_ID == "")) {
    stop("Integration_Batch_ID construction produced missing values")
  }
  meta
}

integration_batch_audit <- function(meta) {
  required <- c(
    "Dataset", "Integration_Batch_ID", "Integration_Batch_Model",
    "Integration_Batch_Source", "Global_Donor_ID", "Specimen_ID",
    "Technical_Batch_ID", "Technical_Batch_Verification_Status",
    "Technical_Batch_Evidence", "Technical_Batch_Use_Status",
    "Technical_Batch_Use_Evidence",
    "Integration_Batch_ID_Dataset",
    "Integration_Batch_ID_Verified_Technical",
    "Technical_Batch_Dataset_Eligible",
    "Technical_Batch_Dataset_Eligibility_Reason"
  )
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L) {
    stop("integration batch audit is missing columns: ", paste(
      missing,
      collapse = ", "
    ))
  }
  groups <- split(
    seq_len(nrow(meta)),
    paste(meta$Dataset, meta$Integration_Batch_ID, sep = "\r")
  )
  do.call(rbind, lapply(groups, function(index) {
    donor <- normalize_missing_metadata(meta$Global_Donor_ID[index])
    specimen <- normalize_missing_metadata(meta$Specimen_ID[index])
    batch_id <- normalize_missing_metadata(
      meta$Technical_Batch_ID[index]
    )
    batch_status <- normalize_missing_metadata(
      meta$Technical_Batch_Verification_Status[index]
    )
    batch_evidence <- normalize_missing_metadata(
      meta$Technical_Batch_Evidence[index]
    )
    batch_use_status <- normalize_missing_metadata(
      meta$Technical_Batch_Use_Status[index]
    )
    batch_use_evidence <- normalize_missing_metadata(
      meta$Technical_Batch_Use_Evidence[index]
    )
    batch_source_column <- normalize_missing_metadata(
      meta$Technical_Batch_Source_Column[index]
    )
    batch_derivation_rule <- normalize_missing_metadata(
      meta$Technical_Batch_Derivation_Rule[index]
    )
    collapse_known <- function(value) {
      value <- sort(unique(stats::na.omit(value)))
      if (length(value) == 0L) NA_character_ else paste(value, collapse = ";")
    }
    data.frame(
      Dataset = unique(meta$Dataset[index]),
      Integration_Batch_ID = unique(meta$Integration_Batch_ID[index]),
      Integration_Batch_Model = unique(
        meta$Integration_Batch_Model[index]
      ),
      Integration_Batch_Source = paste(
        sort(unique(meta$Integration_Batch_Source[index])),
        collapse = ";"
      ),
      Source_Technical_Batch_ID = collapse_known(batch_id),
      Technical_Batch_Verification_Status = collapse_known(batch_status),
      Technical_Batch_Evidence = collapse_known(batch_evidence),
      Technical_Batch_Source_Column = collapse_known(batch_source_column),
      Technical_Batch_Derivation_Rule = collapse_known(
        batch_derivation_rule
      ),
      Technical_Batch_Use_Status = collapse_known(batch_use_status),
      Technical_Batch_Use_Evidence = collapse_known(batch_use_evidence),
      Integration_Batch_ID_Dataset = collapse_known(
        meta$Integration_Batch_ID_Dataset[index]
      ),
      Integration_Batch_ID_Verified_Technical = collapse_known(
        meta$Integration_Batch_ID_Verified_Technical[index]
      ),
      Technical_Batch_Dataset_Eligible = all(
        as.logical(meta$Technical_Batch_Dataset_Eligible[index])
      ),
      Technical_Batch_Dataset_Eligibility_Reason = collapse_known(
        meta$Technical_Batch_Dataset_Eligibility_Reason[index]
      ),
      Cells = length(index),
      Known_Donors = length(unique(stats::na.omit(donor))),
      Known_Specimens = length(unique(stats::na.omit(specimen))),
      stringsAsFactors = FALSE
    )
  }))
}

read_feature_metadata <- function(file) {
  if (!file.exists(file)) {
    stop("Feature metadata is missing: ", file)
  }
  read_brainomics_table(file)
}

feature_crosswalk <- function(
  feature_metadata,
  object_features,
  dataset
) {
  if (!"Original_Feature_ID" %in% names(feature_metadata)) {
    stop(dataset, " feature metadata lacks Original_Feature_ID")
  }
  original <- as.character(feature_metadata$Original_Feature_ID)
  if (!identical(original, object_features)) {
    stop(dataset, " feature metadata differs from matrix feature order")
  }

  canonical <- rep(NA_character_, length(original))
  mapping_source <- rep(NA_character_, length(original))
  for (column in c(
    "gene_name",
    "feature_name",
    "gene_symbol",
    "Original_Feature_ID"
  )) {
    if (!column %in% names(feature_metadata)) {
      next
    }
    values <- trimws(as.character(feature_metadata[[column]]))
    values[
      is.na(values) |
        values == "" |
        toupper(values) == "NA"
    ] <- NA_character_
    selected <- is.na(canonical) & !is.na(values)
    canonical[selected] <- values[selected]
    mapping_source[selected] <- column
  }
  if (anyNA(canonical) || any(canonical == "")) {
    stop(dataset, " contains features without a canonical identifier")
  }
  duplicated_canonical <- unique(canonical[duplicated(canonical)])
  if (length(duplicated_canonical) > 0L) {
    stop(
      dataset,
      " feature mapping is not one-to-one; duplicated identifiers: ",
      paste(head(duplicated_canonical, 10L), collapse = ", ")
    )
  }
  data.frame(
    Dataset = dataset,
    Original_Feature_ID = original,
    Canonical_Feature_ID = canonical,
    Mapping_Source = mapping_source,
    Renamed = original != canonical,
    Canonical_Is_Ensembl = grepl(
      "^ENSG[0-9]+(?:[.][0-9]+)?$",
      canonical
    ),
    stringsAsFactors = FALSE
  )
}

harmonize_reference_features <- function(
  object,
  dataset,
  feature_file
) {
  feature_metadata <- read_feature_metadata(feature_file)
  crosswalk <- feature_crosswalk(
    feature_metadata = feature_metadata,
    object_features = rownames(object),
    dataset = dataset
  )
  counts <- processed_counts(object)
  rownames(counts) <- crosswalk$Canonical_Feature_ID
  harmonized <- create_processed_object(
    counts = counts,
    metadata = object[[]],
    dataset = dataset,
    input_cells = colnames(object),
    input_features = crosswalk$Canonical_Feature_ID
  )
  list(
    object = harmonized,
    crosswalk = crosswalk
  )
}

add_metadata_column <- function(meta, name, value) {
  if (!name %in% names(meta)) {
    meta[[name]] <- value
  }
  meta
}

existing_source_cell_type_spec <- function(dataset) {
  specs <- list(
    AllenM1 = list(
      original = c("cluster_label"),
      label = c("cluster_label"),
      level_1 = c("class_label"),
      level_2 = c("subclass_label"),
      level_3 = c("cluster_label"),
      cluster = c("cell_type_accession_label"),
      source_file = "Allen Brain Map human M1 metadata.csv",
      original_semantics = "original Allen taxonomy cluster label",
      label_semantics = "original Allen taxonomy cluster label"
    ),
    EGAD00001006049 = list(
      original = c("Source_CellType_Original"),
      label = c("Source_CellType_Original"),
      level_1 = c("Source_CellClass"),
      level_2 = c("Source_CellType_Original"),
      cluster = c("Source_Cluster_ID"),
      doublet = c("Source_Doublet_Flag"),
      source_file = "public author HumanFetalBrainPool.h5 tensor",
      original_semantics = "author cluster Class indexed by source Clusters",
      label_semantics = "author cluster Class indexed by source Clusters"
    ),
    EGAS00001006537 = list(
      original = c("Source_CellType_Original"),
      label = c("Source_CellType_Original"),
      level_1 = c("Source_CellType_Original"),
      source_file = "Cameron et al. public Figshare per-nucleus metadata",
      original_semantics = "original author cellIDs annotation",
      label_semantics = "original author cellIDs annotation"
    ),
    GSE144136 = list(
      original = c("Source_CellType_Original"),
      label = c("Source_CellType_Original"),
      level_1 = c("Source_CellType_Original"),
      source_file = "GSE144136 author-annotated CellNames.csv",
      original_semantics = "original author fine cell-type prefix",
      label_semantics = "original author fine cell-type prefix"
    ),
    GSE178175 = list(
      original = c("Source_CellType_Original"),
      label = c("Source_CellType_Original"),
      level_1 = c("Source_CellType_Original"),
      source_file = "GEO GSE178175 original 10x matrices and family SOFT",
      original_semantics = "original study did not release per-nucleus labels",
      label_semantics = "unavailable in the original public study files"
    ),
    GSE202210 = list(
      original = c("Source_CellType_Original"),
      label = c("Source_CellType_Original"),
      level_1 = c("Source_CellType_Original"),
      source_file = paste(
        "GEO GSE202210 original 10x matrices plus the original HCA",
        "project-submission donor metadata"
      ),
      original_semantics = "original study did not release per-nucleus labels",
      label_semantics = "unavailable in the original public study files"
    ),
    GSE103723 = list(
      label = c("cell_type", "CellType", "Cell_type", "CellType_raw"),
      level_1 = c("cell_type", "CellType_raw"),
      level_2 = c("CellType"),
      level_3 = c("Cell_type"),
      cluster = c("priCluster", "cluster1", "seurat_clusters", "cluster")
    ),
    GSE104276 = list(
      original = "cell_types",
      label = "cell_types",
      level_1 = "cell_types",
      level_2 = "Source_Author_Neuronal_Class",
      source_file = "GSE104276_readme_sample_barcode.xlsx: SampleInfo and PFC_ExcitatoryNeurons_rf_c7"
    ),
    GSE186538 = list(
      label = c("original_name", "cell_type", "CellType_raw"),
      level_1 = c("original_name", "CellType_raw"),
      level_2 = c("subclass", "cell_type"),
      cluster = c("seurat_clusters", "cluster")
    ),
    GSE199762 = list(
      label = c("all.exp_type", "cell_type", "CellType_raw"),
      level_1 = c("all.exp_type", "CellType_raw"),
      level_2 = c("cell_type", "CellType"),
      cluster = c("seurat_clusters", "cluster")
    ),
    GSE204683 = list(
      label = c("Cell.type", "CellType_raw"),
      level_1 = c("Cell.type", "CellType_raw"),
      source_file = "GSE204683_barcodes.tsv: exact donor and barcode join (11=EaFet2; 16=EaFet1)"
    ),
    GSE207334 = list(
      label = c("subclass", "cell_type", "CellType_raw"),
      level_1 = c("class", "subclass", "CellType_raw"),
      level_2 = c("subclass", "cell_type"),
      cluster = c("leiden", "seurat_clusters", "cluster")
    ),
    GSE212606 = list(
      label = c("Cell_type", "Cell_Subtype", "CellType_raw"),
      level_1 = c("Cell_type", "CellType_raw"),
      level_2 = c("Cell_Subtype", "cell_subtype"),
      cluster = c("seurat_clusters", "cluster")
    ),
    GSE217511 = list(
      label = c("celltypes", "cell_type", "CellType_raw"),
      level_1 = c("celltypes", "CellType_raw"),
      level_2 = c("cell_type"),
      cluster = c("seurat_clusters", "RNA_snn_res.0.8", "cluster")
    ),
    GSE261983 = list(
      label = c("anno", "subclass", "azimuth", "CellType_raw"),
      level_1 = c("subclass", "CellType_raw"),
      level_2 = c("anno"),
      level_3 = c("azimuth"),
      source_file = paste(
        "brainSCOPE Girgenti-multiome annotated matrices and",
        "Girgenti-multiome_cell_metadata.tsv"
      )
    ),
    GSE296073 = list(
      label = c("cluster1", "cluster2", "CellType_raw"),
      level_1 = c("cluster1", "CellType_raw"),
      level_2 = c("cluster2"),
      cluster = c("leiden", "seurat_clusters", "cluster")
    ),
    GSE67835 = list(
      label = c("cluster", "CellType_raw"),
      level_1 = c("cluster", "CellType_raw")
    ),
    GSE81475 = list(
      original = c("CellType", "CellType_raw"),
      label = c("CellType", "CellType_raw"),
      level_1 = c("CellType", "CellType_raw"),
      cluster = c("cluster1", "cluster", "seurat_clusters"),
      source_file = paste(
        "GSE81475 retained source annotations plus GEO series-matrix",
        "cell-to-GSM metadata"
      )
    ),
    GSE97942 = list(
      label = c("priCluster", "cell_type", "CellType_raw"),
      level_1 = c("priCluster", "CellType_raw"),
      level_2 = c("cell_type"),
      cluster = c("cluster", "seurat_clusters")
    ),
    HYPOMAP = list(
      original = c("C4_named", "celltype_annotation", "CellType_raw"),
      label = c("C4_named", "C3_named", "celltype_annotation", "CellType_raw"),
      level_1 = c("C0_named", "celltype_annotation"),
      level_2 = c("C1_named"),
      level_3 = c("C2_named"),
      cluster = c("C4", "C3", "C2"),
      source_file = paste(
        "public HYPOMAP author H5AD; C0_named-C4_named retain the",
        "published hierarchical neuronal and non-neuronal annotations"
      )
    ),
    Li_et_al_2018 = list(
      original = c("Source_Cell_Type_Fine", "CellType_raw"),
      label = c("Source_Cell_Type_Fine", "CellType_raw"),
      level_1 = c("Source_Cell_Type_Broad", "Source_Cell_Type_Fine"),
      level_2 = c("Source_Cell_Type_Fine"),
      cluster = c("subtype", "ctype"),
      source_file = paste(
        "Li et al. 2018 original fetal and adult PsychENCODE RData;",
        "prenatal Table S3 and adult Table S4"
      )
    ),
    Ma_et_al_2022 = list(
      original = c("subtype", "subclass", "class", "CellType_raw"),
      label = c("subtype", "subclass", "class", "CellType_raw"),
      level_1 = c("class"),
      level_2 = c("subclass"),
      level_3 = c("subtype"),
      cluster = c("seurat_clusters", "cluster")
    ),
    Nowakowski_et_al_2017 = list(
      original = c("WGCNAcluster", "CellType_raw"),
      label = c("WGCNAcluster", "CellType_raw"),
      level_1 = c("WGCNAcluster", "CellType_raw"),
      cluster = c("WGCNAcluster"),
      source_file = "public cortex-dev rawData/meta.tsv"
    ),
    PRJCA015229 = list(
      label = c("Annotation", "BigCellType", "CellType_raw"),
      level_1 = c("BigCellType", "CellType_raw"),
      level_2 = c("Annotation"),
      level_3 = c("SubCellType", "cell_subtype"),
      cluster = c("seurat_clusters", "leiden"),
      doublet = c(
        "Source_Scrublet_Class", "Source_Predicted_Doublet",
        "PredictedDoublets", "db.Scrublet_class",
        "predicted_doublet", "doublet", "Doublet"
      )
    ),
    GSE168408 = list(
      original = c("sub_clust"),
      label = c("sub_clust"),
      level_1 = c("cell_type"),
      level_2 = c("major_clust"),
      level_3 = c("sub_clust"),
      cluster = c("leiden"),
      source_file = paste(
        "Lister Lab RNA-all_full-counts-and-downsampled-CPM.h5ad obs",
        "+ RNA-all_BCs-meta-data.csv"
      ),
      original_semantics =
        "original Herring et al. final subcluster annotation",
      label_semantics =
        "preferred finest original-author annotation"
    ),
    ROSMAP = list(
      label = c("Subcelltype", "Celltype", "CellType_raw"),
      level_1 = c("Celltype", "CellType_raw"),
      level_2 = c("Subcelltype"),
      cluster = c("seurat_clusters", "cluster"),
      doublet = c("predicted_doublet", "doublet")
    ),
    SomaMut = list(
      original = c("new_clusters", "CellType_raw"),
      label = c("new_clusters", "new_clusters2", "new_clusters3"),
      level_1 = c("new_clusters3"),
      level_2 = c("new_clusters2"),
      level_3 = c("new_clusters"),
      cluster = c("seurat_clusters", "cluster")
    )
  )
  spec <- specs[[dataset]]
  if (is.null(spec)) {
    spec <- list(
      label = "CellType_raw",
      level_1 = "CellType_raw"
    )
  }
  if (is.null(spec$source_file)) {
    spec$source_file <- "source metadata retained during dataset preprocessing"
  }
  if (is.null(spec$original)) {
    spec$original <- spec$label
  }
  if (is.null(spec$original_semantics)) {
    spec$original_semantics <-
      "original-study annotation retained without harmonization"
  }
  if (is.null(spec$label_semantics)) {
    spec$label_semantics <-
      "preferred source-study annotation retained without harmonization"
  }
  spec
}

first_reported_metadata_field <- function(meta, candidates) {
  value <- rep(NA_character_, nrow(meta))
  source <- rep(NA_character_, nrow(meta))
  for (candidate in candidates) {
    if (!candidate %in% names(meta)) {
      next
    }
    candidate_value <- normalize_missing_metadata(meta[[candidate]])
    selected <- is.na(value) & !is.na(candidate_value)
    value[selected] <- candidate_value[selected]
    source[selected] <- candidate
  }
  list(name = source, value = value)
}

# Recover author labels by source identity; never substitute atlas predictions.
recover_verified_source_cell_types <- function(meta, dataset) {
  if (!dataset %in% c("GSE204683", "GSE104276")) return(meta)
  cells <- if ("Original_Cell_ID" %in% names(meta)) meta$Original_Cell_ID else meta$Cells
  stopifnot(length(cells) == nrow(meta), !anyNA(cells), !anyDuplicated(cells))
  if (dataset == "GSE204683") {
    source <- data.table::fread(brainomics_data_path("raw", dataset, "GSE204683_barcodes.tsv"))
    donor <- c(`4` = "LaFet1", `8` = "LaFet2", `11` = "EaFet2", `16` = "EaFet1",
      `150656` = "Adult2", `150666` = "Adult1", `4413` = "Inf1", `4422` = "Inf2",
      `5936` = "Adol2", `5977` = "Child2", `6007` = "Adol1", `6032` = "Child1")
    keys <- paste(source[["Donor ID"]], source$Barcode, sep = "::")
    query <- paste(unname(donor[sub("_.*", "", cells)]), sub("^[^_]+_", "", cells), sep = "::")
    at <- match(query, keys)
    stopifnot(!anyDuplicated(keys), !anyNA(at))
    meta$Cell.type <- as.character(source[["Cell type"]][at])
    meta$CellType_raw <- meta$Cell.type
  } else {
    file <- brainomics_data_path("raw", dataset, "rawData", "GSE104276_readme_sample_barcode.xlsx")
    info <- readxl::read_xlsx(file, sheet = "SampleInfo", col_names = FALSE, .name_repair = "minimal")
    exc <- readxl::read_xlsx(file, sheet = "PFC_ExcitatoryNeurons_rf_c7", col_names = FALSE, .name_repair = "minimal")
    ids <- as.character(info[[1]][-1L])
    stopifnot(!anyNA(ids), !anyDuplicated(ids))
    meta$cell_types <- as.character(info[[2]][-1L])[match(cells, ids)]
    meta$CellType_raw <- meta$cell_types
    meta$Source_Cell_Type_Broad_Original <- meta$cell_types
    meta$Source_Author_Neuronal_Class <- ifelse(cells %in% as.character(exc[[1]][-1L]),
      "Excitatory neurons", NA_character_)
  }
  meta
}

add_existing_source_cell_type_metadata <- function(meta, dataset) {
  meta <- recover_verified_source_cell_types(meta, dataset)
  spec <- existing_source_cell_type_spec(dataset)
  selected <- lapply(
    c("original", "label", "level_1", "level_2", "level_3", "cluster"),
    function(level) first_reported_metadata_field(meta, spec[[level]])
  )
  names(selected) <- c(
    "original", "label", "level_1", "level_2", "level_3", "cluster"
  )
  annotation_candidates <- unique(unlist(
    spec[c(
      "original", "label", "level_1", "level_2", "level_3", "cluster"
    )],
    use.names = FALSE
  ))
  available_fields <- annotation_candidates[
    vapply(annotation_candidates, function(field) {
      field %in% names(meta) && any(!is.na(
        normalize_missing_metadata(meta[[field]])
      ))
    }, logical(1))
  ]
  meta$CellType_raw <- selected$label$value
  meta$Source_CellType <- selected$label$value
  meta$source_cell_type_original_label <- selected$original$value
  meta$source_cell_type_label <- selected$label$value
  meta$source_cell_type_level_1 <- selected$level_1$value
  meta$source_cell_type_level_2 <- selected$level_2$value
  meta$source_cell_type_level_3 <- selected$level_3$value
  meta$source_cluster_id <- selected$cluster$value
  meta$source_cell_type_original_label_source_column <-
    selected$original$name
  meta$source_cell_type_original_label_semantics <-
    spec$original_semantics
  meta$source_cell_type_label_semantics <- spec$label_semantics
  meta$source_cell_type_label_source_column <- selected$label$name
  meta$source_cell_type_level_1_source_column <- selected$level_1$name
  meta$source_cell_type_level_2_source_column <- selected$level_2$name
  meta$source_cell_type_level_3_source_column <- selected$level_3$name
  meta$source_cluster_id_source_column <- selected$cluster$name
  meta$source_cell_type_source_column <- if (
    length(available_fields) == 0L
  ) {
    NA_character_
  } else {
    paste(available_fields, collapse = ";")
  }
  meta$source_cell_type_source_file <- spec$source_file
  annotation_available <- !is.na(selected$original$value) |
    !is.na(selected$label$value) |
    !is.na(selected$level_1$value) |
    !is.na(selected$level_2$value) |
    !is.na(selected$level_3$value) |
    !is.na(selected$cluster$value)
  meta$source_cell_type_annotation_status <- ifelse(
    annotation_available,
    "source annotation hierarchy retained",
    "not available"
  )
  meta$source_cell_type_annotation_scope <- ifelse(
    annotation_available,
    paste(
      "available source annotation fields are retained with their exact",
      "source columns; unavailable cells remain explicitly unannotated"
    ),
    "not applicable"
  )
  meta$source_cell_type_mapping_method <- ifelse(
    annotation_available,
    "identity; source labels and cluster IDs retained",
    "no source label available"
  )
  confidence <- first_reported_metadata_field(
    meta,
    spec$confidence
  )$value
  meta$source_cell_type_confidence <- ifelse(
    annotation_available & is.na(confidence),
    "source study assignment; numeric confidence not reported",
    confidence
  )
  meta$source_cell_type_doublet_flag <- first_reported_metadata_field(
    meta,
    spec$doublet
  )$value
  meta
}

existing_source_sex_spec <- function(dataset) {
  candidates <- c(
    "Sex_Source_Raw", "donor_gender", "sex", "Sex", "gender", "Gender"
  )
  if (identical(dataset, "EGAD00001006049")) {
    candidates <- c("Sex_Source_Raw", "Sex_impute", candidates)
  }
  list(candidates = unique(candidates))
}

add_existing_source_sex_metadata <- function(meta, dataset) {
  spec <- existing_source_sex_spec(dataset)
  selected <- first_reported_metadata_field(meta, spec$candidates)
  meta$Sex_Source_Raw <- selected$value
  source_column <- selected$name
  if ("Sex_Source_Column" %in% names(meta)) {
    retained_source_column <- normalize_missing_metadata(
      meta$Sex_Source_Column
    )
    use_retained <- !is.na(retained_source_column)
    source_column[use_retained] <- retained_source_column[use_retained]
  }
  meta$Sex_Source_Column <- source_column
  assignment <- rep("not reported", nrow(meta))
  reported <- !is.na(selected$value)
  assignment[reported] <- "source reported"
  imputed <- reported & !is.na(selected$name) & grepl(
    "imput",
    selected$name,
    ignore.case = TRUE
  )
  assignment[imputed] <- "upstream source imputation; not author inferred"
  if ("Sex_Assignment_Method" %in% names(meta)) {
    retained <- normalize_missing_metadata(meta$Sex_Assignment_Method)
    use_retained <- !is.na(retained)
    assignment[use_retained] <- retained[use_retained]
  }
  meta$Sex_Assignment_Method <- assignment
  meta
}

upgrade_existing_reference_metadata <- function(meta, dataset) {
  stopifnot(is.data.frame(meta), length(dataset) == 1L)
  n <- nrow(meta)
  original_cells <- if ("Cells" %in% names(meta)) {
    as.character(meta$Cells)
  } else {
    rownames(meta)
  }
  source_sample <- if ("Original_Source_Sample_ID" %in% names(meta)) {
    normalize_missing_metadata(meta$Original_Source_Sample_ID)
  } else if ("Original_Sample" %in% names(meta)) {
    normalize_missing_metadata(meta$Original_Sample)
  } else if ("Sample" %in% names(meta)) {
    normalize_missing_metadata(meta$Sample)
  } else {
    rep(NA_character_, n)
  }
  source_record <- if ("Original_Source_Record_ID" %in% names(meta)) {
    normalize_missing_metadata(meta$Original_Source_Record_ID)
  } else if ("Original_Library_ID" %in% names(meta)) {
    normalize_missing_metadata(meta$Original_Library_ID)
  } else if ("Sample_ID" %in% names(meta)) {
    normalize_missing_metadata(meta$Sample_ID)
  } else {
    rep(NA_character_, n)
  }

  donor <- if ("Original_Donor_ID" %in% names(meta)) {
    normalize_missing_metadata(meta$Original_Donor_ID)
  } else {
    source_sample
  }

  specimen <- if ("Original_Specimen_ID" %in% names(meta)) {
    normalize_missing_metadata(meta$Original_Specimen_ID)
  } else {
    source_sample
  }
  library <- rep(NA_character_, n)
  technical_batch <- rep(NA_character_, n)

  donor_rule <- rep(
    "source Sample field interpreted as donor",
    n
  )
  donor_status <- rep("source field interpretation", n)
  specimen_meaning <- rep(
    "source donor-level biological specimen",
    n
  )
  specimen_status <- rep("source field interpretation", n)
  library_meaning <- rep(
    paste(
      "not populated; source record is preserved separately and was not",
      "independently verified as a sequencing library"
    ),
    n
  )
  library_status <- rep("unknown", n)
  library_evidence <- rep(
    paste(
      "no source field has been promoted to Library_ID without",
      "paper- or repository-level verification"
    ),
    n
  )
  library_source_column <- rep(NA_character_, n)
  library_derivation_rule <- rep(NA_character_, n)
  technical_batch_meaning <- rep(
    "not available as an independently verified technical factor",
    n
  )
  technical_batch_status <- rep("not available", n)
  technical_batch_evidence <- rep(
    paste(
      "no paper- or repository-verified library-preparation, sequencing-run,",
      "lane, GEM-well, or chemistry-batch field retained in this object"
    ),
    n
  )
  technical_batch_source_column <- rep(NA_character_, n)
  technical_batch_derivation_rule <- rep(NA_character_, n)
  technical_batch_use_status <- rep("not available", n)
  technical_batch_use_evidence <- rep(
    paste(
      "no eligible technical correction factor was admitted for this",
      "dataset"
    ),
    n
  )

  rebuilt_original_sources <- c(
    "AllenM1", "EGAD00001006049", "EGAS00001006537",
    "GSE144136", "GSE178175", "GSE202210"
  )
  if (dataset %in% rebuilt_original_sources) {
    required_rebuilt <- c(
      "Original_Donor_ID", "Original_Specimen_ID",
      "Original_Library_ID", "Original_Source_Record_ID"
    )
    missing_rebuilt <- setdiff(required_rebuilt, names(meta))
    if (length(missing_rebuilt) > 0L) {
      stop(
        dataset,
        " rebuilt source metadata is missing: ",
        paste(missing_rebuilt, collapse = ", ")
      )
    }
    donor <- normalize_missing_metadata(meta$Original_Donor_ID)
    specimen <- normalize_missing_metadata(meta$Original_Specimen_ID)
    library <- normalize_missing_metadata(meta$Original_Library_ID)
    source_record <- normalize_missing_metadata(meta$Original_Source_Record_ID)
    technical_batch <- if ("Original_Technical_Batch_ID" %in% names(meta)) {
      normalize_missing_metadata(meta$Original_Technical_Batch_ID)
    } else {
      rep(NA_character_, n)
    }
    if (anyNA(donor) || anyNA(specimen) || anyNA(library) ||
      anyNA(source_record)) {
      stop(dataset, " rebuilt biological-unit hierarchy is incomplete")
    }
    donor_rule[] <- "retain the exact author/repository donor identifier"
    donor_status[] <- "verified: original author/repository metadata"
    specimen_meaning[] <- "source anatomical tissue specimen"
    specimen_status[] <- "verified: original author/repository metadata"
    library_meaning[] <- if (dataset == "GSE144136") {
      "donor-level GEO 10x gene-expression library record"
    } else {
      "source 10x gene-expression library or sample record"
    }
    library_status[] <- "verified: original author/repository identifier"
    library_evidence[] <- switch(dataset,
      AllenM1 = paste(
        "Allen Brain Map human M1 metadata identifies the two donor-level",
        "10x gene-expression source records"
      ),
      EGAD00001006049 = paste(
        "HumanFetalBrainPool.h5 SampleID identifies each author-processed",
        "10x library"
      ),
      EGAS00001006537 = paste(
        "Figshare per-nucleus metadata supplies the exact donor-region",
        "sample/library value"
      ),
      GSE144136 = paste(
        "GEO GSE144136 family SOFT maps each donor number to one GSM",
        "10x source record"
      ),
      GSE178175 = paste(
        "GEO GSE178175 supplies two donor-specific 10x matrix GSM records"
      ),
      GSE202210 = paste(
        "GEO GSE202210 supplies 12 donor-specific 10x matrix GSM records"
      )
    )
    library_source_column[] <- "Original_Library_ID"
    library_derivation_rule[] <-
      "retain the exact original-source library/sample identifier"

    has_technical <- !is.na(technical_batch)
    technical_batch_meaning[has_technical] <- if (
      dataset == "EGAD00001006049"
    ) {
      "author-reported 10x chemistry version"
    } else if (dataset == "GSE144136") {
      "author cell-ID B# preparation batch"
    } else if (dataset == "EGAS00001006537") {
      "author sample B# value"
    } else {
      "source library-specific technical record"
    }
    technical_batch_status[has_technical] <-
      "verified: exact original-source field"
    technical_batch_evidence[has_technical] <- switch(dataset,
      EGAD00001006049 = paste(
        "Braun et al. author tensor reports Chemistry for every cell;",
        "the paper used 10x Chromium v2 and v3"
      ),
      GSE144136 = paste(
        "Nagy et al. author cell identifiers encode the B# preparation",
        "batch for every nucleus"
      ),
      EGAS00001006537 = paste(
        "Cameron et al. Figshare metadata reports the B# suffix in every",
        "donor-region sample identifier"
      ),
      "original repository source record retained without reinterpretation"
    )
    technical_batch_source_column[has_technical] <-
      "Original_Technical_Batch_ID"
    technical_batch_derivation_rule[has_technical] <-
      "retain the exact source field"
    analysis_include <- if ("Analysis_Include" %in% names(meta)) {
      as.logical(meta$Analysis_Include)
    } else {
      rep(TRUE, n)
    }
    analysis_has_technical <- has_technical & analysis_include
    batch_donor_count <- if (any(analysis_has_technical)) {
      vapply(
        split(
          donor[analysis_has_technical],
          technical_batch[analysis_has_technical]
        ),
        function(value) length(unique(value)),
        integer(1)
      )
    } else {
      integer()
    }
    eligible_batch_design <- dataset %in% c(
      "EGAD00001006049", "GSE144136"
    ) && length(batch_donor_count) >= 2L && all(batch_donor_count >= 2L)
    technical_batch_use_status[has_technical] <- if (
      eligible_batch_design
    ) {
      "eligible for technical-batch sensitivity model"
    } else {
      "record only: ineligible for technical-batch correction"
    }
    technical_batch_use_evidence[has_technical] <- if (
      eligible_batch_design
    ) {
      paste(
        "the verified field has at least two levels among analysis-eligible",
        "cells and every retained level spans multiple donors"
      )
    } else {
      paste(
        "at least one source technical level is donor-specific among",
        "analysis-eligible cells or the field has fewer than two usable",
        "analysis levels"
      )
    }
  }

  specimen_from_source_record <- dataset %in% c(
    "HYPOMAP", "GSE199762", "GSE103723",
    "Li_et_al_2018", "GSE207334",
    "GSE97942"
  )
  if (specimen_from_source_record) {
    specimen <- source_record
    specimen_meaning[] <-
      "source sample/specimen record nested within donor"
  }

  if (dataset == "GSE103723") {
    source_pattern <- "^.+_([0-9]+W[FM])_B[0-9]+$"
    valid_source_label <- !is.na(source_sample) & grepl(
      source_pattern,
      source_sample
    )
    if (!all(valid_source_label)) {
      stop("GSE103723 source region-donor-batch labels are malformed")
    }
    donor <- sub(source_pattern, "\\1", source_sample)
    specimen <- sub("_B[0-9]+$", "", source_sample)
    if (!identical(sort(unique(donor)), c("22WF", "23WF", "23WM"))) {
      stop("GSE103723 must retain the three published embryo donors")
    }
    donor_rule[] <- paste(
      "extract the age-sex embryo key 22WF, 23WF, or 23WM from the",
      "source region_age-sex_B# label; the paper reports one 22-week",
      "embryo and two 23-week embryos"
    )
    donor_status[] <- "curated from publication-defined embryo identity"
    specimen_meaning[] <-
      "embryo-by-dissected-brain-region biological specimen"
    specimen_status[] <- "curated from source region and embryo label"

    library <- source_record
    if (anyNA(library) || length(unique(library)) != 51L) {
      stop("GSE103723 must retain all 51 GEO library records")
    }
    library_meaning[] <-
      "STRT-seq library pooling up to 96 individually barcoded cells"
    library_status[] <- "verified: GEO library accession"
    library_evidence[] <- paste(
      "Fan et al. 2018, DOI 10.1038/s41422-018-0053-3, reports that",
      "cDNA from all 96 differently barcoded cells was pooled for one",
      "library construction; GEO GSE103723 supplies 51 corresponding GSM",
      "records for the retained 4,664 cells"
    )
    library_source_column[] <- "Sample_ID"
    library_derivation_rule[] <-
      "retain the exact GEO GSM accession stored in Sample_ID"

    technical_batch <- library
    technical_batch_meaning[] <-
      "source STRT-seq pooled-cell library"
    technical_batch_status[] <-
      "verified: publication library protocol and GEO records"
    technical_batch_evidence[] <- library_evidence
    technical_batch_source_column[] <- "Sample_ID"
    technical_batch_derivation_rule[] <- library_derivation_rule
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each of the 51 library levels belongs to one embryo-region",
      "specimen and one donor; the cross-donor design gate fails"
    )
  }

  if (dataset == "GSE217511") {
    donor <- sub("[CG]$", "", source_sample)
    specimen <- source_sample
    donor_rule[] <- paste(
      "paired region suffix C/G removed from source sample;",
      "supported by GEO sample pairing"
    )
    donor_status[] <- "curated from explicit paired-region source labels"
    specimen_meaning[] <-
      "donor-region specimen encoded by source C/G suffix"
    specimen_status[] <- "source reported"
    library <- source_record
    library_meaning[] <- paste(
      "10x 3-prime gene-expression library for one donor-region specimen;",
      "identified by its one-to-one GEO GSM accession"
    )
    library_status[] <- "verified: source library identifier"
    library_evidence[] <- paste(
      "Ramos et al. 2022, DOI 10.1038/s41467-022-34975-2, states that",
      "each region and sample was barcoded and sequenced as a separate",
      "library; GEO GSE217511 provides one GSM matrix record for each",
      "source sample code"
    )
    library_source_column[] <- "Sample_ID"
    library_derivation_rule[] <- paste(
      "retain the GEO geo_accession stored in Sample_ID after verifying",
      "a one-to-one mapping with the source donor-region sample code"
    )
    mapping <- unique(data.frame(
      library = library,
      specimen = specimen,
      stringsAsFactors = FALSE
    ))
    if (anyNA(mapping) || anyDuplicated(mapping$library) ||
      anyDuplicated(mapping$specimen)) {
      stop("GSE217511 GSM-to-specimen library mapping is not one-to-one")
    }

    technical_batch <- library
    technical_batch_meaning[] <-
      "donor-region-specific 10x gene-expression library"
    technical_batch_status[] <- paste(
      "verified: GSE217511 GEO library records and publication methods"
    )
    technical_batch_evidence[] <- library_evidence
    technical_batch_source_column[] <- "Sample_ID"
    technical_batch_derivation_rule[] <- library_derivation_rule
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each library contains one donor-region specimen, so library and",
      "biological unit are confounded and the cross-donor design gate fails"
    )
  }

  if (dataset == "GSE204683" &&
    any(!is.na(source_record))) {
    donor <- source_record
    donor_rule[] <-
      "source numeric donor identifier retained from Sample_ID"
    donor_status[] <- "source reported"
    library <- source_record
    library_meaning[] <- paste(
      "one 10x Multiome gene-expression library generated for the source",
      "brain sample; identified by the retained numeric source sample ID"
    )
    library_status[] <- "verified: one library per source brain sample"
    library_evidence[] <- paste(
      "Trevino et al. 2023, DOI 10.1126/sciadv.adg3754, states that one",
      "library was generated for each brain sample; GSE204683 retains the",
      "numeric donor/sample identifier in Sample_ID"
    )
    library_source_column[] <- "Sample_ID"
    library_derivation_rule[] <- paste(
      "use the retained numeric source sample ID as the library key because",
      "the publication documents exactly one library per brain sample"
    )
    mapping <- unique(data.frame(
      library = library,
      specimen = specimen,
      stringsAsFactors = FALSE
    ))
    if (anyNA(mapping) || anyDuplicated(mapping$library) ||
      anyDuplicated(mapping$specimen)) {
      stop("GSE204683 sample-to-library mapping is not one-to-one")
    }

    technical_batch <- library
    technical_batch_meaning[] <-
      "donor/sample-specific 10x Multiome gene-expression library"
    technical_batch_status[] <-
      "verified: one library per source brain sample"
    technical_batch_evidence[] <- library_evidence
    technical_batch_source_column[] <- "Sample_ID"
    technical_batch_derivation_rule[] <- library_derivation_rule
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each verified library contains one donor/sample and is confounded",
      "with the biological unit; the CELLxGENE Batch field is preserved",
      "as a candidate but is not promoted because its technical semantics",
      "are not defined by the paper or public repository metadata"
    )
  }

  if (dataset == "GSE168408") {
    if (!"batch" %in% names(meta)) {
      stop("GSE168408 source metadata is missing the audited batch field")
    }
    library <- normalize_missing_metadata(meta$batch)
    if (anyNA(library)) {
      stop("GSE168408 source library identifiers are incomplete")
    }
    mapping <- unique(data.frame(
      library = library,
      donor = donor,
      stringsAsFactors = FALSE
    ))
    library_donors <- tapply(
      mapping$donor,
      mapping$library,
      function(value) length(unique(stats::na.omit(value)))
    )
    if (any(library_donors != 1L)) {
      stop("GSE168408 source library maps to more than one donor")
    }
    library_meaning[] <- paste(
      "source snRNA-seq library identifier combining sample, age label,",
      "and 10x chemistry version"
    )
    library_status[] <- "verified: source RNA library identifier"
    library_evidence[] <- paste(
      "the author H5AD batch field identifies 27 RNA libraries by donor,",
      "age and chemistry and maps one-to-one to GEO snRNA-seq records"
    )
    library_source_column[] <- "batch"
    library_derivation_rule[] <- "retain the exact author H5AD batch field"

    if (!"Herring_Technical_Batch_Source" %in% names(meta)) {
      stop("GSE168408 lacks the audited platform-chemistry field")
    }
    technical_batch <- normalize_missing_metadata(
      meta$Herring_Technical_Batch_Source
    )
    if (anyNA(technical_batch)) {
      stop("GSE168408 platform-chemistry identifiers are incomplete")
    }
    technical_batch_donors <- tapply(
      donor,
      technical_batch,
      function(value) length(unique(stats::na.omit(value)))
    )
    if (length(technical_batch_donors) != 3L ||
      any(technical_batch_donors < 2L)) {
      stop("GSE168408 technical batches do not span expected donors")
    }
    technical_batch_meaning[] <-
      "cross-donor sequencing-platform and 10x chemistry combination"
    technical_batch_status[] <- paste(
      "verified: author chemistry field and GEO platform records"
    )
    technical_batch_evidence[] <- paste(
      "the author H5AD retains v2/v3 chemistry and GEO assigns each RNA",
      "library to NextSeq 550 or NovaSeq 6000"
    )
    technical_batch_source_column[] <-
      "Herring_Technical_Batch_Source"
    technical_batch_derivation_rule[] <- paste(
      "Source_GEO_Platform_ID + '_' + author chem"
    )
    technical_batch_use_status[] <-
      "eligible for technical-batch sensitivity model"
    technical_batch_use_evidence[] <- paste(
      "all three platform-chemistry combinations span at least two donors;",
      "Dataset remains the primary study-level batch"
    )
  }

  if (dataset == "GSE207334") {
    if (!"Source_GEO_RNA_Record_ID" %in% names(meta)) {
      stop(
        "GSE207334 source metadata is missing the audited RNA GSM field"
      )
    }
    library <- normalize_missing_metadata(meta$Source_GEO_RNA_Record_ID)
    if (anyNA(library)) {
      stop("GSE207334 RNA GSM library identifiers are incomplete")
    }
    mapping <- unique(data.frame(
      library = library,
      donor = donor,
      stringsAsFactors = FALSE
    ))
    if (anyDuplicated(mapping$library) || anyDuplicated(mapping$donor)) {
      stop("GSE207334 RNA library-to-donor mapping is not one-to-one")
    }
    library_meaning[] <-
      "GEO RNA library from one donor-level 10x Multiome experiment"
    library_status[] <- "verified: GEO RNA library accession"
    library_evidence[] <- paste(
      "GEO GSE207334 records GSM6284664-GSM6284668 as the five",
      "transcriptomic libraries generated with Chromium Single Cell",
      "Multiome ATAC plus Gene Expression and sequenced on NovaSeq 6000"
    )
    library_source_column[] <- "Source_GEO_RNA_Record_ID"
    library_derivation_rule[] <-
      "retain the exact donor-to-RNA-GSM mapping from GEO GSE207334"

    technical_batch <- library
    technical_batch_meaning[] <-
      "donor-specific 10x Multiome gene-expression library"
    technical_batch_status[] <-
      "verified: GSE207334 GEO RNA library records"
    technical_batch_evidence[] <- library_evidence
    technical_batch_source_column[] <- "Source_GEO_RNA_Record_ID"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each of the five RNA library levels contains exactly one donor,",
      "so the prespecified cross-donor design gate fails"
    )
  }

  if (dataset == "Ma_et_al_2022") {
    if (!"tech_rep" %in% names(meta)) {
      stop("Ma et al. 2022 metadata is missing the source tech_rep field")
    }
    library <- normalize_missing_metadata(meta$tech_rep)
    if (anyNA(library) ||
      any(sub("_[1-4]$", "", library) != donor)) {
      stop("Ma et al. 2022 technical replicates do not map to donors")
    }
    library_meaning[] <-
      "source 10x snRNA-seq technical-replicate library"
    library_status[] <- "verified: source technical-replicate identifier"
    library_evidence[] <- paste(
      "DOI 10.1126/science.abo7257 states that the human species group",
      "contains four donors with four technical replicates per donor;",
      "the source processed object retains all 16 exact tech_rep values"
    )
    library_source_column[] <- "tech_rep"
    library_derivation_rule[] <- "retain the exact source tech_rep field"

    technical_batch <- library
    technical_batch_meaning[] <-
      "donor-specific 10x snRNA-seq technical replicate"
    technical_batch_status[] <- paste(
      "verified: Science 2022 methods and source tech_rep metadata"
    )
    technical_batch_evidence[] <- library_evidence
    technical_batch_source_column[] <- "tech_rep"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each technical-replicate level contains one donor even though each",
      "donor has four replicates; the prespecified cross-donor gate fails"
    )
  }

  if (dataset == "GSE296073") {
    region <- if ("BrainRegion" %in% names(meta)) {
      normalize_missing_metadata(meta$BrainRegion)
    } else if ("Brain_Region" %in% names(meta)) {
      normalize_missing_metadata(meta$Brain_Region)
    } else {
      rep(NA_character_, n)
    }
    specimen <- ifelse(
      is.na(source_sample) | is.na(region),
      NA_character_,
      paste(source_sample, region, sep = "::")
    )
    library <- source_record
    specimen_meaning[] <-
      "source donor-by-brain-region biological specimen"
    library_meaning[] <- "source libraryID sequencing library"
    library_status[] <- "verified source library identifier"
    library_evidence[] <- paste(
      "GSE296073 source processed object field libraryID contains explicit",
      "library records (for example h2023001whgw22_all1);",
      "DOI 10.1038/s41586-025-09362-8 and GEO GSE296073"
    )
    library_source_column[] <- "libraryID"
    library_derivation_rule[] <-
      "retain the exact source processed-object libraryID field"

    technical_batch <- library
    technical_batch_meaning[] <-
      "donor- and tissue/cell-sort-specific source sequencing library"
    technical_batch_status[] <-
      "verified: source libraryID and GSE296073 publication/repository"
    technical_batch_evidence[] <- library_evidence
    technical_batch_source_column[] <- "libraryID"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each source library is nested within one donor and a biological",
      "tissue or cell-sort stratum, so the cross-donor design gate fails"
    )
  }

  if (dataset == "GSE261983") {
    specimen <- source_sample
    library <- source_record
    specimen_meaning[] <-
      "donor-level prefrontal-cortex tissue specimen"
    specimen_status[] <- "source reported"
    library_meaning[] <- paste(
      "GEO RNA library/source record for one donor-level paired-multiome",
      "prefrontal-cortex specimen"
    )
    library_status[] <- "verified: GEO RNA library/source record"
    library_evidence[] <- paste(
      "GEO GSE261983 deposits one GSM*_RT*_RNA record per RNA library;",
      "the processing input retains this exact value in",
      "Source_GEO_Record_ID and Sample_ID; DOI 10.1126/science.adi5199"
    )
    library_source_column[] <- "Source_GEO_Record_ID"
    library_derivation_rule[] <-
      "retain the exact Source_GEO_Record_ID RNA record as Library_ID"
    mapping <- unique(data.frame(
      library = library,
      specimen = specimen,
      stringsAsFactors = FALSE
    ))
    if (anyNA(mapping) || anyDuplicated(mapping$library) ||
      anyDuplicated(mapping$specimen)) {
      stop("GSE261983 RNA-library-to-specimen mapping is not one-to-one")
    }

    technical_batch <- library
    technical_batch_meaning[] <-
      "donor-specific paired-multiome RNA library/source record"
    technical_batch_status[] <- "verified: GSE261983 GEO RNA record"
    technical_batch_evidence[] <- library_evidence
    technical_batch_source_column[] <- "Source_GEO_Record_ID"
    technical_batch_derivation_rule[] <- library_derivation_rule
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each RNA library contains one donor-level specimen, so library and",
      "donor are confounded and the cross-donor design gate fails"
    )
  }

  if (dataset == "GSE97942") {
    required_source_fields <- c(
      "Source_GEO_Record_ID", "Source_Experiment_ID", "Source_Donor_ID"
    )
    missing_source_fields <- setdiff(required_source_fields, names(meta))
    if (length(missing_source_fields) > 0L) {
      stop(
        "GSE97942 is missing audited GEO fields: ",
        paste(missing_source_fields, collapse = ", ")
      )
    }
    donor <- normalize_missing_metadata(meta$Source_Donor_ID)
    region <- if ("Region" %in% names(meta)) {
      normalize_missing_metadata(meta$Region)
    } else {
      normalize_missing_metadata(meta$Brain_Region)
    }
    specimen <- ifelse(
      is.na(donor) | is.na(region),
      NA_character_,
      paste(donor, region, sep = "::")
    )
    library <- normalize_missing_metadata(meta$Source_GEO_Record_ID)
    technical_batch <- normalize_missing_metadata(meta$Source_Experiment_ID)
    if (anyNA(donor) || anyNA(specimen) || anyNA(library) ||
      anyNA(technical_batch) || length(unique(donor)) != 6L ||
      length(unique(library)) != 46L ||
      length(unique(technical_batch)) != 20L) {
      stop("GSE97942 donor-library-experiment hierarchy is incomplete")
    }
    experiment_donors <- vapply(
      split(donor, technical_batch),
      function(value) length(unique(value)),
      integer(1)
    )
    if (any(experiment_donors != 1L)) {
      stop("GSE97942 experiment spans more than one donor")
    }
    donor_rule[] <- paste(
      "retain GEO patient identifiers where reported; two anonymous",
      "donors are grouped from the exact GEO age, sex and region records",
      "under the paper's six-donor design"
    )
    donor_status[] <-
      "verified or curated anonymous group from GEO sample records"
    specimen_meaning[] <- "donor-by-brain-region biological specimen"
    specimen_status[] <- "curated from GEO donor and region"
    library_meaning[] <-
      "sample-indexed snDrop-seq library deposited as one GEO GSM record"
    library_status[] <- "verified: GEO RNA library accession"
    library_evidence[] <- paste(
      "Lake et al. 2018, DOI 10.1038/nbt.4038, reports 46 libraries",
      "prepared from 20 experiments across six individuals; GSE97930",
      "MINiML maps every retained library title to its GSM record"
    )
    library_source_column[] <- "Source_GEO_Record_ID"
    library_derivation_rule[] <- paste(
      "map the D7_ orig.ident library title to the exact GSE97930 GSM",
      "using the frozen 46-record GEO MINiML crosswalk"
    )
    technical_batch_meaning[] <-
      "source snDrop-seq experiment in which one to six libraries were split"
    technical_batch_status[] <-
      "verified: publication methods and GSE97930 experiment metadata"
    technical_batch_evidence[] <- paste(
      "the paper reports 20 experiments split into 46 indexed libraries;",
      "GSE97930 supplies the experiment date/batch for each GSM record"
    )
    technical_batch_source_column[] <- "Source_Experiment_ID"
    technical_batch_derivation_rule[] <-
      "retain the exact GSE97930 experiment characteristic"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each of the 20 experiments contains exactly one donor and region;",
      "the cross-donor design gate fails"
    )
  }

  if (dataset == "PRJCA015229") {
    required_prjca <- c(
      "Source_RNA_Library_ID", "Source_Scrublet_Score",
      "Source_Scrublet_Class", "BigCellType", "Annotation"
    )
    missing_prjca <- setdiff(required_prjca, names(meta))
    if (length(missing_prjca) > 0L) {
      stop(
        "PRJCA015229 is missing audited source fields: ",
        paste(missing_prjca, collapse = ", ")
      )
    }
    library <- normalize_missing_metadata(meta$Source_RNA_Library_ID)
    if (anyNA(library) || length(unique(library)) != 6L) {
      stop("PRJCA015229 must retain six human RNA libraries")
    }
    mapping <- unique(data.frame(
      library = library,
      donor = donor,
      stringsAsFactors = FALSE
    ))
    if (anyNA(mapping) || anyDuplicated(mapping$library) ||
      anyDuplicated(mapping$donor)) {
      stop("PRJCA015229 RNA-library-to-donor mapping is not one-to-one")
    }
    specimen <- donor
    specimen_meaning[] <-
      "donor-level anterior cingulate cortex biological specimen"
    specimen_status[] <- "source reported"
    library_meaning[] <- paste(
      "donor-specific RNA library from either standalone snRNA-seq or the",
      "gene-expression arm of single-nucleus Multiome"
    )
    library_status[] <- "verified: source assay and donor identifier"
    library_evidence[] <- paste(
      "Yuan et al. 2024, DOI 10.1016/j.xgen.2024.100703, reports four",
      "human single-nucleus Multiome samples and two standalone human",
      "snRNA-seq samples; PRJCA015229 metadata retains Sample and Tech"
    )
    library_source_column[] <- "Source_RNA_Library_ID"
    library_derivation_rule[] <-
      "concatenate exact source Sample and Tech fields with ::"

    technical_batch <- library
    technical_batch_meaning[] <-
      "donor-specific RNA assay library"
    technical_batch_status[] <-
      "verified: Cell Genomics 2024 methods and PRJCA015229 metadata"
    technical_batch_evidence[] <- paste(
      library_evidence,
      "RNA libraries were sequenced on Illumina NovaSeq 6000."
    )
    technical_batch_source_column[] <- "Source_RNA_Library_ID"
    technical_batch_derivation_rule[] <- library_derivation_rule
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each of the six RNA libraries contains exactly one donor and assay",
      "type is partly confounded with donor; the cross-donor gate fails"
    )
  }

  if (dataset == "GSE212606") {
    if (!"EasySci_PCR_Group_ID" %in% names(meta)) {
      stop(
        "GSE212606 source metadata is missing the audited ",
        "EasySci_PCR_Group_ID field"
      )
    }
    pcr_group <- normalize_missing_metadata(meta$EasySci_PCR_Group_ID)
    expected_group <- sub("[.][^.]+$", "", original_cells)
    if (anyNA(pcr_group) || !identical(pcr_group, expected_group)) {
      stop(
        "GSE212606 EasySci PCR groups do not match the Cell_ID prefix"
      )
    }
    donor_count <- vapply(
      split(donor, pcr_group),
      function(value) length(unique(stats::na.omit(value))),
      integer(1)
    )
    if (length(donor_count) < 2L || any(donor_count < 2L)) {
      stop(
        "GSE212606 retained PCR groups must each span multiple donors"
      )
    }

    library <- rep("GSM6657986", n)
    library_meaning[] <- paste(
      "GEO aggregate human EasySci-RNA library/source record containing",
      "the non-demultiplexed hippocampus and SMTG experiment"
    )
    library_status[] <- "verified: GEO EasySci-RNA library record"
    library_evidence[] <- paste(
      "GEO GSM6657986 identifies the human experiment as an EasySci-RNA",
      "library and reports Illumina NovaSeq 6000; the record does not",
      "claim donor-level libraries"
    )
    library_source_column[] <- "GEO accession constant"
    library_derivation_rule[] <- paste(
      "assign constant GSM6657986 only as the aggregate GEO library/source",
      "record; do not reinterpret donor IDs as libraries"
    )

    technical_batch <- pcr_group
    technical_batch_meaning[] <- paste(
      "final indexed EasySci PCR group (P7-indexed final PCR well)",
      "retained as the FASTQ-demultiplexing group"
    )
    technical_batch_status[] <- paste(
      "verified: GEO GSM6657986 protocol, Supplementary Protocol 1, and",
      "official EasySci pipeline Zenodo 8395492"
    )
    technical_batch_evidence[] <- paste(
      "GEO GSM6657986 documents redistribution to final PCR wells, addition",
      "of an indexed P7 primer per well, and pooling of final PCR products;",
      "the official EasySci pipeline defines each demultiplexed FASTQ",
      "sample ID as a PCR-group ID; the exact group is the source Cell_ID",
      "prefix before the ligation-plus-RT barcode"
    )
    technical_batch_source_column[] <- "EasySci_PCR_Group_ID"
    technical_batch_derivation_rule[] <- paste(
      "remove the final dot and 20-base ligation-plus-RT barcode from",
      "Cell_ID; the resulting exact source field is the PCR-group ID"
    )
    technical_batch_use_status[] <-
      "eligible for technical-batch sensitivity model"
    technical_batch_use_evidence[] <- paste(
      "the retained cohort has complete PCR-group metadata and every level",
      "spans multiple control donors; the full 30,801-cell source audit",
      "confirms 120 levels and all six donors in every level"
    )
  }

  if (dataset == "SomaMut") {
    required_somamut <- c(
      "orig.ident", "case", "batch", "new_clusters3",
      "new_clusters2", "new_clusters", "predicted.id",
      "prediction.score.max"
    )
    missing_somamut <- setdiff(required_somamut, names(meta))
    if (length(missing_somamut) > 0L) {
      stop(
        "SomaMut source metadata is missing audited fields: ",
        paste(missing_somamut, collapse = ", ")
      )
    }
    donor <- normalize_missing_metadata(meta$orig.ident)
    specimen <- donor
    library <- normalize_missing_metadata(meta$case)
    if (anyNA(donor) || anyNA(library)) {
      stop("SomaMut donor or library identifiers are missing")
    }
    library_donor <- sub("_.*$", "", library)
    if (!identical(library_donor, donor)) {
      stop("SomaMut library identifiers do not map exactly to orig.ident donors")
    }
    donor_rule[] <- "retain the exact source processed-object orig.ident field"
    donor_status[] <- "source reported"
    specimen_meaning[] <-
      "donor-level prefrontal-cortex tissue specimen"
    specimen_status[] <- "source reported"
    library_meaning[] <- paste(
      "source snRNA-seq library or sorted-library replicate identified by",
      "the exact public pfc.clean.rds case value"
    )
    library_status[] <- "verified: source library identifier"
    library_evidence[] <- paste(
      "Jeffries et al. 2025 Supplementary Table 3 labels all 42 retained",
      "case values as libraries and reports their post-QC cell counts;",
      "the counts sum to the 367,317 nuclei in DOI",
      "10.1038/s41586-025-09435-8"
    )
    library_source_column[] <- "case"
    library_derivation_rule[] <-
      "retain the exact public pfc.clean.rds case field"

    technical_batch <- normalize_missing_metadata(meta$batch)
    if (anyNA(technical_batch)) {
      stop("SomaMut exact source batch is missing")
    }
    library_batch <- unique(data.frame(
      library = library,
      technical_batch = technical_batch,
      stringsAsFactors = FALSE
    ))
    if (anyDuplicated(library_batch$library)) {
      stop("SomaMut source library maps to multiple preparation batches")
    }
    technical_batch_meaning[] <- paste(
      "snRNA-seq preparation batch reported for each retained library in",
      "Supplementary Table 3"
    )
    technical_batch_status[] <- paste(
      "verified: Jeffries et al. Nature 2025 methods and Supplementary",
      "Table 3"
    )
    technical_batch_evidence[] <- paste(
      "DOI 10.1038/s41586-025-09435-8 states that samples were prepared",
      "in batches of up to six donors with mixed ages and sexes;",
      "Supplementary Table 3 gives a batch value for all 42 retained",
      "libraries and all 367,317 nuclei"
    )
    technical_batch_source_column[] <- "batch"
    technical_batch_derivation_rule[] <-
      "retain the exact public pfc.clean.rds batch value"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "although the source batch mapping is complete, at least one retained",
      "batch level contains only one donor; therefore not every batch",
      "crosses biological units and the prespecified design gate fails"
    )
  }

  if (dataset == "GSE199762") {
    nonmultiplexed_donors <- c("H69", "H71", "H31", "H37")
    multiplexed_donors <- c("H48", "H39", "H46", "H29", "H33")
    known_design_donors <- c(nonmultiplexed_donors, multiplexed_donors)
    if (any(!is.na(donor) & !donor %in% known_design_donors)) {
      stop("GSE199762 contains a donor outside the published run design")
    }
    technical_batch <- ifelse(
      donor %in% nonmultiplexed_donors,
      "experiment_1_nonmultiplexed",
      ifelse(
        donor %in% multiplexed_donors,
        "experiment_2_multiplexed",
        NA_character_
      )
    )
    technical_batch_meaning[] <- paste(
      "published nonmultiplexed versus multiplexed experimental run;",
      "derived from the explicitly listed source donors"
    )
    technical_batch_status[] <- paste(
      "verified: Nascimento et al. Nature 2024 methods and GEO GSE199762"
    )
    technical_batch_evidence[] <- paste(
      "DOI 10.1038/s41586-023-06981-x describes one nonmultiplexed run",
      "for H69/H71/H31/H37 and a separate multiplexed run for",
      "H48/H39/H46/H29/H33"
    )
    technical_batch_source_column[] <- "Sample"
    technical_batch_derivation_rule[] <- paste(
      "H69,H71,H31,H37=experiment_1_nonmultiplexed;",
      "H48,H39,H46,H29,H33=experiment_2_multiplexed"
    )
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "the two experimental runs have intentionally different donor-age",
      "and region composition, including concentration of early-childhood",
      "donors in the multiplexed run; correction would risk removing biology"
    )
  }

  if (dataset == "GSE104276") {
    required_gse104276 <- c(
      "Original_Source_Sample_ID", "Source_GEO_Record_ID",
      "Source_Single_Cell_Library_ID", "Source_Published_Cell_Count",
      "Source_Retained_Cell_Count"
    )
    if (!all(required_gse104276 %in% names(meta))) {
      stop("GSE104276 audited donor, GEO and library fields are required")
    }
    donor <- normalize_missing_metadata(meta$Original_Source_Sample_ID)
    source_gsm <- normalize_missing_metadata(meta$Source_GEO_Record_ID)
    single_cell_library <- normalize_missing_metadata(
      meta$Source_Single_Cell_Library_ID
    )
    if (
      anyNA(donor) || anyNA(source_gsm) || anyNA(single_cell_library) ||
        length(unique(donor)) != 12L ||
        length(unique(source_gsm)) != 39L ||
        length(unique(single_cell_library)) != 2392L ||
        !all(meta$Source_Published_Cell_Count == 2394L) ||
        !all(meta$Source_Retained_Cell_Count == 2392L)
    ) {
      stop("GSE104276 audited 12-donor/39-GSM/2,392-library design differs")
    }
    gsm_donors <- vapply(
      split(donor, source_gsm),
      function(value) length(unique(value)),
      integer(1)
    )
    if (any(gsm_donors != 1L)) {
      stop("GSE104276 GEO records must remain nested within one donor")
    }
    specimen <- donor
    donor_rule[] <- paste(
      "derive the publication biological sample from the",
      "GW##_PFC# prefix of the source cell identifier"
    )
    donor_status[] <- paste(
      "verified: publication/GEO design has 12 biological PFC samples"
    )
    specimen_meaning[] <- "publication-defined donor PFC specimen"
    specimen_status[] <- "verified: GEO sample design"

    library <- single_cell_library
    library_meaning[] <-
      "individual modified-Smart-seq2 single-cell RNA library"
    library_status[] <- "verified: source single-cell library identifier"
    library_evidence[] <- paste(
      "Zhong et al. 2018, DOI 10.1038/nature25980, used a modified",
      "Smart-seq2 single-cell protocol with UMI barcodes; GEO GSE104276",
      "retains the exact cell/library labels and 39 pooled GSM records"
    )
    library_source_column[] <- "Source_Single_Cell_Library_ID"
    library_derivation_rule[] <- "retain the exact source cell/library ID"

    technical_batch <- source_gsm
    technical_batch_meaning[] <-
      "GEO sequencing record pooling multiple single-cell libraries"
    technical_batch_status[] <-
      "verified: 39 GEO sequencing records reconstructed from raw archive"
    technical_batch_evidence[] <- paste(
      library_evidence,
      "The raw GEO archive maps the 2,392 retained cells to 39 GSM records."
    )
    technical_batch_source_column[] <- "Source_GEO_Record_ID"
    technical_batch_derivation_rule[] <- paste(
      "map exact cell headers from the 39 GEO raw files; resolve ten",
      "published boundary labels using the source barcode workbook"
    )
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "every GSM record contains exactly one of the 12 donor specimens;",
      "the cross-donor design gate fails. The source paper reports 2,394",
      "cells whereas the retained matrix contains 2,392, recorded as",
      "source-to-input attrition rather than silently changing the count."
    )
  }
  if (dataset == "GSE67835") {
    required_gse67835 <- c(
      "Original_Donor_ID", "Original_Specimen_ID",
      "Source_GEO_Record_ID", "Source_Single_Cell_Library_ID",
      "Source_Capture_Device_ID", "Source_Instrument_Model",
      "Source_Full_Series_Cell_Count", "Source_Retained_Cell_Count"
    )
    missing_gse67835 <- setdiff(required_gse67835, names(meta))
    if (length(missing_gse67835) > 0L) {
      stop(
        "GSE67835 audited donor, library and platform fields are missing: ",
        paste(missing_gse67835, collapse = ", ")
      )
    }
    donor <- normalize_missing_metadata(meta$Original_Donor_ID)
    specimen <- normalize_missing_metadata(meta$Original_Specimen_ID)
    library <- normalize_missing_metadata(
      meta$Source_Single_Cell_Library_ID
    )
    source_record <- normalize_missing_metadata(meta$Source_GEO_Record_ID)
    capture_device <- normalize_missing_metadata(
      meta$Source_Capture_Device_ID
    )
    technical_batch <- normalize_missing_metadata(
      meta$Source_Instrument_Model
    )
    if (
      n != 466L || anyNA(donor) || anyNA(specimen) || anyNA(library) ||
        anyNA(source_record) || anyNA(capture_device) ||
        anyNA(technical_batch) || length(unique(donor)) != 12L ||
        length(unique(specimen)) != 14L ||
        length(unique(library)) != 466L ||
        length(unique(source_record)) != 466L ||
        length(unique(capture_device)) != 16L ||
        length(unique(technical_batch)) != 2L ||
        !all(meta$Source_Full_Series_Cell_Count == 466L) ||
        !all(meta$Source_Retained_Cell_Count == 466L)
    ) {
      stop("GSE67835 audited 466-cell/12-donor/16-chip design differs")
    }
    donor_rule[] <- paste(
      "retain GEO experiment_sample_name as the publication biological",
      "donor key (8 adult AB and 4 prenatal FB specimens)"
    )
    donor_status[] <- "verified: cell-level GEO donor/sample field"
    specimen_meaning[] <- "donor-by-source-tissue biological specimen"
    specimen_status[] <- "verified: GEO donor and tissue fields"
    library_meaning[] <- "individual Fluidigm C1 single-cell RNA library"
    library_status[] <- "verified: one GEO GSM record per single cell"
    library_evidence[] <- paste(
      "Darmanis et al. 2015, DOI 10.1073/pnas.1507125112, reports 466",
      "retained single cells; GEO GSE67835 exposes one GSM count file",
      "per published cell"
    )
    library_source_column[] <- "Source_Single_Cell_Library_ID"
    library_derivation_rule[] <- "retain the exact cell-level GEO GSM"
    technical_batch_meaning[] <-
      "GEO-reported sequencing instrument model"
    technical_batch_status[] <-
      "verified: cell-level GEO instrument_model field"
    technical_batch_evidence[] <- paste(
      "GSE67835 reports Illumina MiSeq or Illumina NextSeq 500 for every",
      "cell; each platform level contains cells from six source donors"
    )
    technical_batch_source_column[] <- "Source_Instrument_Model"
    technical_batch_derivation_rule[] <-
      "retain the exact GEO instrument_model value"
    technical_batch_use_status[] <-
      "eligible for technical-batch sensitivity model"
    technical_batch_use_evidence[] <- paste(
      "both verified platform levels span six donors and multiple ages;",
      "platform is used only in the prespecified technical sensitivity",
      "model, while the 16 donor-confounded C1 chips remain audit fields"
    )
  }

  if (dataset == "GSE81475") {
    required_gse81475 <- c(
      "Original_Donor_ID", "Original_Specimen_ID",
      "Source_GEO_Record_ID", "Source_Single_Cell_Library_ID",
      "Source_Capture_Device_ID", "Source_Sequencing_Lane",
      "Source_Full_Series_Cell_Count", "Source_Retained_Cell_Count"
    )
    missing_gse81475 <- setdiff(required_gse81475, names(meta))
    if (length(missing_gse81475) > 0L) {
      stop(
        "GSE81475 audited metadata fields are missing: ",
        paste(missing_gse81475, collapse = ", ")
      )
    }
    donor <- normalize_missing_metadata(meta$Original_Donor_ID)
    specimen <- normalize_missing_metadata(meta$Original_Specimen_ID)
    library <- normalize_missing_metadata(
      meta$Source_Single_Cell_Library_ID
    )
    source_sample <- specimen
    source_record <- normalize_missing_metadata(meta$Source_GEO_Record_ID)
    technical_batch <- normalize_missing_metadata(
      meta$Source_Capture_Device_ID
    )
    sequencing_lane <- normalize_missing_metadata(
      meta$Source_Sequencing_Lane
    )
    if (
      n != 476L || anyNA(donor) || anyNA(specimen) || anyNA(library) ||
        anyNA(source_record) || anyNA(technical_batch) ||
        anyNA(sequencing_lane) || length(unique(donor)) != 3L ||
        length(unique(library)) != 476L ||
        length(unique(source_record)) != 476L ||
        length(unique(technical_batch)) != 6L ||
        !all(meta$Source_Full_Series_Cell_Count == 1608L) ||
        !all(meta$Source_Retained_Cell_Count == 476L)
    ) {
      stop("GSE81475 audited 476-cell/3-donor/6-chip design differs")
    }
    chip_donors <- vapply(
      split(donor, technical_batch),
      function(value) length(unique(value)),
      integer(1)
    )
    if (any(chip_donors != 1L)) {
      stop("GSE81475 Fluidigm C1 chips must remain donor-confounded")
    }
    donor_rule[] <- "retain the donor field from GEO GSE81475 characteristics"
    donor_status[] <- "verified: GEO donor characteristic"
    specimen_meaning[] <- "donor dorsolateral-prefrontal-cortex specimen"
    specimen_status[] <- "verified: GEO donor, age and fetal DFC design"
    library_meaning[] <- paste(
      "individual Fluidigm C1/Nextera XT single-cell RNA library with",
      "one-to-one GEO GSM record"
    )
    library_status[] <- "verified: exact GEO single-cell library title"
    library_evidence[] <- paste(
      "GEO GSE81475 exposes 1,608 cell-level records; the retained",
      "476-cell non-infected fetal DFC subset maps one-to-one to 476 GSM",
      "records and exact single-cell library titles"
    )
    library_source_column[] <- "Source_Single_Cell_Library_ID"
    library_derivation_rule[] <- "retain the exact GEO Sample_title"
    technical_batch_meaning[] <- "Fluidigm C1 capture chip"
    technical_batch_status[] <- paste(
      "verified: GEO characteristic c1_chip; sequencing lane retained",
      "separately in Source_Sequencing_Lane"
    )
    technical_batch_evidence[] <- paste(
      "GEO GSE81475 provides donor, age and c1_chip for every GSM;",
      "all 476 retained cell titles also encode their sequencing lane"
    )
    technical_batch_source_column[] <- "Source_Capture_Device_ID"
    technical_batch_derivation_rule[] <-
      "retain GEO characteristic_c1_chip without concatenating covariates"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "each of the six C1 chips contains exactly one donor; chip and",
      "biological unit are confounded. Lane labels are not promoted because",
      "they are embedded cell-level records without a globally unique run ID"
    )
  }

  if (dataset == "Li_et_al_2018") {
    required_li <- c(
      "Source_Li_Arm", "Original_Donor_ID", "Original_Specimen_ID",
      "Original_Library_ID", "Original_Technical_Batch_ID",
      "Source_Cell_Type_Fine"
    )
    missing_li <- setdiff(required_li, names(meta))
    if (length(missing_li) > 0L) {
      stop(
        "Li et al. 2018 audited metadata fields are missing: ",
        paste(missing_li, collapse = ", ")
      )
    }
    arm <- normalize_missing_metadata(meta$Source_Li_Arm)
    donor <- normalize_missing_metadata(meta$Original_Donor_ID)
    specimen <- normalize_missing_metadata(meta$Original_Specimen_ID)
    library <- normalize_missing_metadata(meta$Original_Library_ID)
    technical_batch <- normalize_missing_metadata(
      meta$Original_Technical_Batch_ID
    )
    source_sample <- specimen
    source_record <- normalize_missing_metadata(
      meta$Original_Source_Record_ID
    )
    adult <- arm == "adult_snRNA"
    fetal <- arm == "fetal_scRNA"
    if (
      n != 17855L || anyNA(arm) || anyNA(donor) || anyNA(specimen) ||
        anyNA(library) || anyNA(technical_batch) || anyNA(source_record) ||
        sum(adult) != 17093L || sum(fetal) != 762L ||
        length(unique(donor[adult])) != 3L ||
        length(unique(donor[fetal])) != 8L ||
        length(unique(library[adult])) != 3L ||
        length(unique(library[fetal])) != 762L ||
        length(unique(technical_batch[adult])) != 3L ||
        length(unique(technical_batch[fetal])) != 12L
    ) {
      stop("Li et al. audited 17,093-adult/762-fetal design differs")
    }
    batch_donors <- vapply(
      split(donor, paste(arm, technical_batch, sep = "::")),
      function(value) length(unique(value)),
      integer(1)
    )
    if (any(batch_donors != 1L)) {
      stop("Li et al. source library/chip levels must remain donor-confounded")
    }
    donor_rule[] <- paste(
      "retain Braincode/Donor from original prenatal meta2 and remove the",
      "DFC suffix from original adult orig.ident"
    )
    donor_status[] <- "verified: original PsychENCODE RData and Tables S3-S4"
    specimen_meaning[adult] <- "adult donor DFC tissue specimen"
    specimen_meaning[fetal] <- "prenatal donor-region tissue specimen"
    specimen_status[] <- "verified: original paper source metadata"
    library_meaning[adult] <- "donor-level 10x single-nucleus RNA library"
    library_meaning[fetal] <-
      "individual Fluidigm C1 single-cell RNA library"
    library_status[] <- "verified: original paper source files"
    library_evidence[adult] <- paste(
      "Li et al. 2018 Table S4 and adult RData retain three adult DFC",
      "sample/library codes for 17,093 nuclei"
    )
    library_evidence[fetal] <- paste(
      "Li et al. 2018 prenatal RData retains 762 exact cell/library IDs;",
      "every ID maps one-to-one to a GSE81475 GSM record"
    )
    library_source_column[adult] <- "Source_10x_Library_ID"
    library_source_column[fetal] <- "Source_Single_Cell_Library_ID"
    library_derivation_rule[adult] <-
      "map adult Braincode to the exact Table S4 Cellcode"
    library_derivation_rule[fetal] <-
      "retain the exact prenatal source cell/library identifier"
    technical_batch_meaning[adult] <-
      "adult donor-specific 10x snRNA-seq library"
    technical_batch_meaning[fetal] <- "prenatal Fluidigm C1 capture chip"
    technical_batch_status[] <-
      "verified: original PsychENCODE RData and paper sample tables"
    technical_batch_evidence[adult] <- library_evidence[adult]
    technical_batch_evidence[fetal] <- paste(
      "the prenatal source identifier records 12 Fluidigm C1 capture chips;",
      "each chip maps to exactly one donor"
    )
    technical_batch_source_column[adult] <- "Source_10x_Library_ID"
    technical_batch_source_column[fetal] <- "Source_Capture_Device_ID"
    technical_batch_derivation_rule[adult] <-
      "retain the adult Table S4 Cellcode"
    technical_batch_derivation_rule[fetal] <-
      "extract the C1-## capture-device prefix from the source cell ID"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "every adult 10x library and prenatal C1 chip contains one donor;",
      "all verified technical levels fail the cross-donor design gate"
    )
  }

  if (dataset == "GSE186538") {
    if (!"batch" %in% names(meta)) {
      stop("GSE186538 source metadata is missing the audited batch field")
    }
    region <- if ("Region" %in% names(meta)) {
      normalize_missing_metadata(meta$Region)
    } else {
      normalize_missing_metadata(meta$subregion)
    }
    specimen <- ifelse(
      is.na(donor) | is.na(region),
      NA_character_,
      paste(donor, region, sep = "::")
    )
    specimen_meaning[] <-
      "donor-by-hippocampal-or-entorhinal-subregion specimen"
    specimen_status[] <- "curated from source donor and subregion"

    library <- normalize_missing_metadata(meta$batch)
    if (anyNA(library) || length(unique(library)) != 25L) {
      stop("GSE186538 must retain all 25 source technical replicates")
    }
    library_donors <- vapply(
      split(donor, library),
      function(value) length(unique(stats::na.omit(value))),
      integer(1)
    )
    if (any(library_donors != 1L)) {
      stop("GSE186538 technical replicate spans more than one donor")
    }
    library_meaning[] <-
      "source sample-indexed 10x snRNA-seq technical-replicate library"
    library_status[] <- "verified: source technical-replicate identifier"
    library_evidence[] <- paste(
      "Franjic et al. 2022, DOI 10.1016/j.neuron.2021.10.036, reports",
      "one to eight technical replicates per human donor and construction",
      "of sample-indexed libraries; the source batch field retains all 25",
      "replicate identifiers"
    )
    library_source_column[] <- "batch"
    library_derivation_rule[] <- "retain the exact source batch field"

    technical_batch <- library
    technical_batch_meaning[] <-
      "donor-specific 10x snRNA-seq technical replicate"
    technical_batch_status[] <-
      "verified: Neuron 2022 methods and source batch metadata"
    technical_batch_evidence[] <- paste(
      library_evidence,
      "The paper further reports mixing uniquely indexed samples across",
      "several HiSeq 4000 lanes to avoid lane bias, so lane is not",
      "reconstructed as a separate cell-level factor."
    )
    technical_batch_source_column[] <- "batch"
    technical_batch_derivation_rule[] <-
      "retain the exact source donor_replicate batch field"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "all 25 technical-replicate levels contain exactly one donor;",
      "the cross-donor design gate fails"
    )
  }
  if (dataset == "HYPOMAP") {
    if (!all(c(
      "Source_Study_Dataset", "Sample_ID", "Donor_ID"
    ) %in% names(meta))) {
      stop("HYPOMAP source study, sample and donor fields are required")
    }
    source_study <- normalize_missing_metadata(meta$Source_Study_Dataset)
    if (
      anyNA(source_study) ||
        !identical(unique(source_study), "Tadross") ||
        n != 311964L
    ) {
      stop(
        "HYPOMAP must retain only the 311,964 newly generated Tadross nuclei"
      )
    }
    donor <- normalize_missing_metadata(meta$Donor_ID)
    source_sample_id <- normalize_missing_metadata(meta$Sample_ID)
    if (
      anyNA(donor) || anyNA(source_sample_id) ||
        length(unique(donor)) != 8L ||
        length(unique(source_sample_id)) != 58L
    ) {
      stop("HYPOMAP Tadross donor or source-sample counts differ")
    }
    sample_donors <- vapply(
      split(donor, source_sample_id),
      function(value) length(unique(value)),
      integer(1)
    )
    if (any(sample_donors != 1L)) {
      stop("HYPOMAP source Sample_ID must remain nested within donor")
    }

    specimen <- source_sample_id
    specimen_meaning[] <-
      "source hypothalamic tissue/nuclei sample nested within donor"
    specimen_status[] <- "verified: publication-defined sample identifier"
    library[] <- NA_character_
    library_meaning[] <- paste(
      "not populated; the paper reports 58 samples and 10x library",
      "generation but does not provide a separate one-to-one library ID"
    )
    library_status[] <- "not independently identifiable"
    library_evidence[] <- paste(
      "Tadross et al. 2025, DOI 10.1038/s41586-024-08504-8, reports",
      "58 newly generated samples, 10x Chromium 3' v3.1 library",
      "generation and NovaSeq 6000 sequencing; Sample_ID is retained as",
      "a sample, not relabelled as an independently verified library"
    )

    technical_batch <- source_sample_id
    technical_batch_meaning[] <- paste(
      "source sample identifier used by the original study as the scVI",
      "batch_key; it is a biological-sample covariate, not a verified",
      "cross-donor sequencing batch"
    )
    technical_batch_status[] <-
      "verified: publication-used sample covariate, not sequencing batch"
    technical_batch_evidence[] <- paste(
      library_evidence,
      "The paper explicitly used Sample_ID as scVI batch_key."
    )
    technical_batch_source_column[] <- "Sample_ID"
    technical_batch_derivation_rule[] <-
      "retain the exact source Sample_ID without concatenation"
    technical_batch_use_status[] <-
      "record only: ineligible for technical-batch correction"
    technical_batch_use_evidence[] <- paste(
      "all 58 Sample_ID levels belong to exactly one donor, so the",
      "cross-donor design gate fails; the Siletti arm is excluded as",
      "overlapping previously published input"
    )
  }
  if (dataset == "Nowakowski_et_al_2017") {
    required_nowakowski <- c(
      "Original_Donor_ID", "Original_Specimen_ID",
      "Source_Single_Cell_Library_ID", "Original_Source_Record_ID",
      "Source_Cell_ID_Prefix", "Source_Full_Cell_Count",
      "Source_Retained_Cell_Count", "WGCNAcluster",
      "Original_File_Type", "Expression_Units",
      "Genome_Transcriptome_Build", "Raw_Integer_Counts_Available"
    )
    missing_nowakowski <- setdiff(required_nowakowski, names(meta))
    if (length(missing_nowakowski) > 0L) {
      stop(
        "Nowakowski audited source fields are missing: ",
        paste(missing_nowakowski, collapse = ", ")
      )
    }
    donor <- normalize_missing_metadata(meta$Original_Donor_ID)
    specimen <- normalize_missing_metadata(meta$Original_Specimen_ID)
    library <- normalize_missing_metadata(
      meta$Source_Single_Cell_Library_ID
    )
    source_sample <- donor
    source_record <- normalize_missing_metadata(
      meta$Original_Source_Record_ID
    )
    source_prefix <- normalize_missing_metadata(
      meta$Source_Cell_ID_Prefix
    )
    technical_batch[] <- NA_character_
    if (
      n != 4261L || anyNA(donor) || anyNA(specimen) || anyNA(library) ||
        anyNA(source_record) || length(unique(donor)) != 48L ||
        length(unique(specimen)) != 92L ||
        length(unique(library)) != 4261L ||
        length(unique(source_record)) != 4261L ||
        sum(!is.na(source_prefix)) != 4232L ||
        length(unique(stats::na.omit(source_prefix))) != 86L ||
        !all(meta$Source_Full_Cell_Count == 4261L) ||
        !all(meta$Source_Retained_Cell_Count == 4261L) ||
        !all(meta$Expression_Units == "TPM") ||
        !all(meta$Genome_Transcriptome_Build == "hg38") ||
        any(meta$Raw_Integer_Counts_Available)
    ) {
      stop("Nowakowski audited 4,261-cell/48-donor/92-specimen design differs")
    }
    prefix_donors <- vapply(
      split(
        donor[!is.na(source_prefix)],
        source_prefix[!is.na(source_prefix)]
      ),
      function(value) length(unique(value)),
      integer(1)
    )
    if (!identical(
      as.integer(table(prefix_donors)),
      c(71L, 11L, 4L)
    )) {
      stop("Nowakowski source cell-ID prefix design differs from audit")
    }
    donor_rule[] <-
      "retain the exact public cortex-dev Name field (Sample1-Sample48)"
    donor_status[] <- paste(
      "verified: public source metadata identifies patient/sample of origin"
    )
    specimen_meaning[] <- paste(
      "donor-by-source RegionName/Laminae/Area biological specimen"
    )
    specimen_status[] <-
      "verified: exact public source anatomical fields"
    library_meaning[] <-
      "individual Fluidigm C1 single-cell RNA expression library"
    library_status[] <-
      "verified: exact source cell/library identifier"
    library_evidence[] <- paste(
      "Nowakowski et al. 2017, DOI 10.1126/science.aap8809, reports",
      "Fluidigm C1 single-cell mRNA sequencing; the public cortex-dev",
      "meta.tsv and TPM matrix contain the same 4,261 exact cell/library IDs"
    )
    library_source_column[] <- "Source_Single_Cell_Library_ID"
    library_derivation_rule[] <- "retain the exact public source Cell value"
    technical_batch_meaning[] <- paste(
      "not populated; a plate-like source cell-ID prefix is preserved in",
      "Source_Cell_ID_Prefix but is not promoted to a technical batch"
    )
    technical_batch_status[] <-
      "not available as an independently verified technical factor"
    technical_batch_evidence[] <- paste(
      library_evidence,
      "The 86 exact prefixes cover 4,232 cells; 71 occur in one donor, 11",
      "in two donors and four in three donors."
    )
    technical_batch_source_column[] <- NA_character_
    technical_batch_derivation_rule[] <- NA_character_
    technical_batch_use_status[] <-
      "not available"
    technical_batch_use_evidence[] <- paste(
      "the field is incomplete for 29 cells, its technical meaning is not",
      "independently verified and 71 of 86 levels are donor-specific; it",
      "fails the prespecified gate and Dataset remains the correction key"
    )
  }

  if (dataset == "ROSMAP") {
    required_rosmap <- c("Batch", "Celltype", "Subcelltype")
    missing_rosmap <- setdiff(required_rosmap, names(meta))
    if (length(missing_rosmap) > 0L) {
      stop(
        "ROSMAP source metadata is missing audited fields: ",
        paste(missing_rosmap, collapse = ", ")
      )
    }
    rosmap_rna_diagnosis <- first_reported_metadata_field(
      meta, c("Diagnosis_RNA_Group_raw", "ADdiag3types")
    )
    rosmap_pathology <- first_reported_metadata_field(
      meta, c("Diagnosis_Pathology_raw", "Pathologic_diagnosis_of_AD")
    )
    if (anyNA(rosmap_rna_diagnosis$value) ||
      anyNA(rosmap_pathology$value) ||
      any(rosmap_rna_diagnosis$value != "nonAD") ||
      any(rosmap_pathology$value != "no")) {
      stop("ROSMAP retained cells must satisfy both non-AD criteria")
    }
    meta$Diagnosis_raw <- paste0(
      "ADdiag3types=", rosmap_rna_diagnosis$value,
      "; Pathologic_diagnosis_of_AD=", rosmap_pathology$value
    )
    technical_batch <- normalize_missing_metadata(meta$Batch)
    if (anyNA(technical_batch)) {
      stop("ROSMAP source Batch is missing for retained cells")
    }
    technical_batch_meaning[] <- paste(
      "source RNA processing/sequencing batch retained in meta.tsv;",
      "each retained non-AD batch spans multiple donors"
    )
    technical_batch_status[] <- paste(
      "verified: Mathys et al. Cell 2023 methods and source meta.tsv"
    )
    technical_batch_evidence[] <- paste(
      "DOI 10.1016/j.cell.2023.08.039 reports individual-batch",
      "alignment and Cell Ranger aggregation; the four retained non-AD",
      "Batch values each span 8-13 retained donors in RNA/meta.tsv"
    )
    technical_batch_source_column[] <- "Batch"
    technical_batch_derivation_rule[] <-
      "retain the exact RNA/meta.tsv Batch value without concatenation"
    batch_donor_count <- vapply(
      split(donor, technical_batch),
      function(value) length(unique(value)),
      integer(1)
    )
    if (any(batch_donor_count < 2L)) {
      stop("ROSMAP each retained Batch must span multiple donors")
    }
    donor_rule[] <- "retain the exact RNA/meta.tsv Individual identifier"
    donor_status[] <- "source reported"
    specimen <- donor
    specimen_meaning[] <-
      "donor-level postmortem prefrontal-cortex biological specimen"
    specimen_status[] <- "verified: source Individual and paper brain region"
    library[] <- NA_character_
    library_meaning[] <- paste(
      "not populated; the released metadata identify donor and aggregate",
      "processing Batch but do not expose a separate RNA library ID"
    )
    library_status[] <- "not independently identifiable"
    library_evidence[] <- paste(
      "Mathys et al. 2023 reports 10x Chromium 3-prime v3 libraries and",
      "Cell Ranger aggregation; public RNA/meta.tsv contains Individual",
      "and Batch but no distinct library identifier"
    )
    technical_batch_use_status[] <-
      "eligible for technical-batch sensitivity model"
    technical_batch_use_evidence[] <- paste(
      "the complete retained Batch field has four levels and every level",
      "spans 8-13 non-AD donors; it is used only in the prespecified",
      "technical sensitivity model"
    )
  }

  assay_type <- first_reported_metadata_field(
    meta,
    c("Assay_Type", "Sequence")
  )$value
  sequencing_platform <- first_reported_metadata_field(
    meta,
    c(
      "Sequencing_Platform", "sequencing_platform", "Platform",
      "platform", "Instrument", "instrument_model"
    )
  )$value
  library_chemistry <- first_reported_metadata_field(
    meta,
    c("Library_Chemistry", "Chemistry", "chemistry")
  )$value
  if (dataset == "GSE212606") {
    meta$Technology <- rep("EasySci-RNA", n)
    meta$Sequence <- rep("snRNA-seq", n)
    meta$Assay_Type <- assay_type <- rep("snRNA-seq", n)
    meta$Sequencing_Platform <- sequencing_platform <- rep(
      "Illumina NovaSeq 6000",
      n
    )
    meta$Library_Chemistry <- library_chemistry <- rep(
      "EasySci-RNA combinatorial indexing workflow",
      n
    )
  }
  meta <- add_metadata_column(
    meta,
    "Original_Cell_ID",
    original_cells
  )
  meta <- add_metadata_column(
    meta,
    "Original_Donor_ID",
    donor
  )
  meta <- add_metadata_column(
    meta,
    "Original_Sample_ID",
    specimen
  )
  meta <- add_metadata_column(
    meta,
    "Original_Specimen_ID",
    specimen
  )
  meta <- add_metadata_column(
    meta,
    "Original_Library_ID",
    library
  )
  meta <- add_metadata_column(
    meta,
    "Original_Source_Sample_ID",
    source_sample
  )
  meta <- add_metadata_column(
    meta,
    "Original_Source_Record_ID",
    source_record
  )
  meta <- add_metadata_column(
    meta,
    "Original_Technical_Batch_ID",
    technical_batch
  )
  meta <- add_metadata_column(
    meta,
    "Donor_ID",
    canonical_dataset_id(dataset, "donor", donor)
  )
  meta <- add_metadata_column(
    meta,
    "Specimen_ID",
    canonical_dataset_id(dataset, "specimen", specimen)
  )
  meta <- add_metadata_column(
    meta,
    "Library_ID",
    canonical_dataset_id(dataset, "library", library)
  )
  meta <- add_metadata_column(
    meta,
    "Source_Record_ID",
    canonical_dataset_id(dataset, "source_record", source_record)
  )
  meta <- add_metadata_column(
    meta,
    "Technical_Batch_ID",
    canonical_dataset_id(dataset, "technical_batch", technical_batch)
  )
  meta <- add_metadata_column(
    meta,
    "Global_Donor_ID",
    canonical_dataset_id(dataset, "donor", donor)
  )
  meta <- add_metadata_column(
    meta,
    "Donor_ID_Biological_Meaning",
    ifelse(
      is.na(donor),
      "unknown; no auditable donor identifier retained",
      "human tissue donor"
    )
  )
  meta <- add_metadata_column(
    meta,
    "Specimen_ID_Biological_Meaning",
    specimen_meaning
  )
  meta <- add_metadata_column(
    meta,
    "Library_ID_Biological_Meaning",
    library_meaning
  )
  meta <- add_metadata_column(
    meta,
    "Source_Record_ID_Biological_Meaning",
    ifelse(
      is.na(source_record),
      "unknown",
      "source sample, accession, batch, project, or cell-level record"
    )
  )
  meta <- add_metadata_column(
    meta,
    "Technical_Batch_ID_Biological_Meaning",
    technical_batch_meaning
  )
  meta <- add_metadata_column(
    meta,
    "Donor_ID_Verification_Status",
    donor_status
  )
  meta <- add_metadata_column(
    meta,
    "Specimen_ID_Verification_Status",
    specimen_status
  )
  meta <- add_metadata_column(
    meta,
    "Library_ID_Verification_Status",
    library_status
  )
  meta <- add_metadata_column(
    meta,
    "Library_ID_Evidence",
    library_evidence
  )
  meta <- add_metadata_column(
    meta,
    "Library_ID_Source_Column",
    library_source_column
  )
  meta <- add_metadata_column(
    meta,
    "Library_ID_Derivation_Rule",
    library_derivation_rule
  )
  meta <- add_metadata_column(
    meta,
    "Technical_Batch_Verification_Status",
    technical_batch_status
  )
  meta <- add_metadata_column(
    meta,
    "Technical_Batch_Evidence",
    technical_batch_evidence
  )
  meta <- add_metadata_column(
    meta,
    "Technical_Batch_Source_Column",
    technical_batch_source_column
  )
  meta <- add_metadata_column(
    meta,
    "Technical_Batch_Derivation_Rule",
    technical_batch_derivation_rule
  )
  meta <- add_metadata_column(
    meta,
    "Technical_Batch_Use_Status",
    technical_batch_use_status
  )
  meta <- add_metadata_column(
    meta,
    "Technical_Batch_Use_Evidence",
    technical_batch_use_evidence
  )
  meta <- add_metadata_column(
    meta,
    "Donor_ID_Derivation_Rule",
    donor_rule
  )
  meta <- add_metadata_column(
    meta,
    "Species",
    rep("Homo sapiens", n)
  )
  meta <- add_metadata_column(meta, "Assay_Type", assay_type)
  meta <- add_metadata_column(
    meta,
    "Sequencing_Platform",
    sequencing_platform
  )
  meta <- add_metadata_column(
    meta,
    "Library_Chemistry",
    library_chemistry
  )
  diagnosis <- first_reported_metadata_field(
    meta,
    c(
      "Diagnosis_raw", "Diagnosis", "diagnosis", "Disorder",
      "Condition", "condition", "disease", "pathology"
    )
  )
  meta$Diagnosis_raw <- diagnosis$value
  meta$Diagnosis_Source_Column <- diagnosis$name
  meta$Diagnosis_Verification_Status <- ifelse(
    is.na(diagnosis$value),
    paste(
      "diagnosis was not retained in this processed object;",
      "source-level inclusion audit required"
    ),
    "source metadata field retained without reinterpretation"
  )
  meta <- add_metadata_column(
    meta,
    "Analysis_Role",
    rep("reference", n)
  )
  meta <- add_metadata_column(
    meta,
    "Analysis_Include",
    rep(TRUE, n)
  )
  meta <- add_metadata_column(
    meta,
    "Exclusion_Reason",
    rep(NA_character_, n)
  )
  meta <- add_metadata_column(
    meta,
    "Duplicate_Group",
    rep(NA_character_, n)
  )
  meta <- add_metadata_column(
    meta,
    "Duplicate_Evidence",
    rep(NA_character_, n)
  )
  meta <- add_existing_source_cell_type_metadata(meta, dataset)
  meta <- add_existing_source_sex_metadata(meta, dataset)
  meta <- apply_existing_dataset_age_rules(meta, dataset)
  meta$Dataset <- dataset
  meta$Reference_Component <- "existing_reference"
  meta$Sample <- meta$Specimen_ID
  meta$Sample_ID <- meta$Specimen_ID
  meta
}

canonicalize_existing_reference_metadata_frame <- function(
  metadata,
  cells,
  dataset,
  donor_crosswalk_file = NULL,
  source_record = NULL
) {
  if (!is.data.frame(metadata) ||
    length(cells) != nrow(metadata) ||
    !identical(rownames(metadata), cells)) {
    stop("existing reference metadata and cell order differ for ", dataset)
  }
  meta <- upgrade_existing_reference_metadata(metadata, dataset)
  if (!is.null(donor_crosswalk_file)) {
    meta <- apply_donor_crosswalk(meta, donor_crosswalk_file)
  }
  if (!is.null(source_record)) {
    meta <- add_source_publication_metadata(meta, source_record)
  }
  meta <- add_age_schema(meta)
  meta <- add_metadata_schema(meta)
  meta$Analysis_Include <- as.logical(meta$Analysis_Include)
  explicit_region_exclusion <- !is.na(meta$BrainRegion) &
    meta$BrainRegion == "Choroid"
  meta$Analysis_Include[explicit_region_exclusion] <- FALSE
  meta$Analysis_Role[explicit_region_exclusion] <-
    "reference_excluded"
  meta$Exclusion_Reason[explicit_region_exclusion] <-
    "choroid plexus excluded from brain parenchyma reference"
  egad_head_exclusion <- dataset == "EGAD00001006049" &
    !is.na(meta$BrainRegion) & meta$BrainRegion == "Head" &
    meta$Analysis_Include
  meta$Analysis_Include[egad_head_exclusion] <- FALSE
  meta$Analysis_Role[egad_head_exclusion] <- "reference_excluded"
  meta$Exclusion_Reason[egad_head_exclusion] <- paste(
    "source Head compartment excluded from the brain reference;",
    "predominantly cranial fibroblast, neural-crest and placodal cells"
  )
  if (dataset == "PRJCA015229") {
    prjca_doublet <- normalize_missing_metadata(
      meta$Source_Scrublet_Class
    )
    if (anyNA(prjca_doublet) ||
      !all(prjca_doublet %in% c("singlet", "doublet"))) {
      stop("PRJCA015229 source Scrublet class is incomplete or invalid")
    }
    doublet_exclusion <- prjca_doublet == "doublet"
    meta$Analysis_Include[doublet_exclusion] <- FALSE
    meta$Analysis_Role[doublet_exclusion] <- "reference_excluded"
    meta$Exclusion_Reason[doublet_exclusion] <-
      "source Scrublet doublet excluded from primary reference analysis"
  }
  meta$Cells <- cells
  rownames(meta) <- cells
  validate_dataset_metadata(
    meta,
    matrix_cells = cells,
    expected_cells = length(cells)
  )
  meta
}

canonicalize_existing_reference_metadata <- function(
  object,
  dataset,
  donor_crosswalk_file = NULL,
  source_record = NULL
) {
  if (!inherits(object, "Seurat")) {
    stop("existing reference input is not a Seurat object")
  }
  canonicalize_existing_reference_metadata_frame(
    metadata = object[[]],
    cells = colnames(object),
    dataset = dataset,
    donor_crosswalk_file = donor_crosswalk_file,
    source_record = source_record
  )
}

read_existing_canonical_metadata <- function(
  file,
  object,
  dataset
) {
  if (!file.exists(file)) {
    stop("canonical metadata sidecar is missing: ", file)
  }
  meta <- read_brainomics_table(file)
  if (any(as.character(meta$Dataset) != dataset) ||
    anyDuplicated(as.character(meta$Cells))) {
    stop(
      "canonical metadata sidecar differs from the source object for ",
      dataset
    )
  }
  object_index <- match(colnames(object), as.character(meta$Cells))
  if (anyNA(object_index)) {
    stop(
      "canonical metadata does not cover every object cell for ",
      dataset
    )
  }
  meta <- meta[object_index, , drop = FALSE]
  if (!identical(as.character(meta$Cells), colnames(object))) {
    stop("canonical metadata cell ordering failed for ", dataset)
  }
  meta$Analysis_Include <- as.logical(meta$Analysis_Include)
  rownames(meta) <- meta$Cells
  validate_dataset_metadata(
    meta,
    matrix_cells = colnames(object),
    expected_cells = ncol(object)
  )
  meta
}

qualify_object_cells <- function(object, dataset) {
  meta <- object[[]]
  if (!identical(rownames(meta), colnames(object))) {
    stop("object metadata row order does not match cells for ", dataset)
  }
  original_cells <- if ("Original_Cell_ID" %in% names(meta)) {
    normalize_missing_metadata(meta$Original_Cell_ID)
  } else {
    colnames(object)
  }
  if (anyNA(original_cells) || anyDuplicated(original_cells)) {
    stop("source cell IDs are missing or duplicated within ", dataset)
  }
  qualified <- paste(dataset, original_cells, sep = "::")
  object <- SeuratObject::RenameCells(object, new.names = qualified)
  meta <- object[[]]
  meta$Cells <- qualified
  meta$Integration_Cell_ID <- qualified
  meta$Original_Cell_ID <- original_cells
  meta$Dataset <- dataset
  rownames(meta) <- qualified
  object@meta.data <- meta
  if (!identical(colnames(object), meta$Cells)) {
    stop("qualified object and metadata cell IDs differ for ", dataset)
  }
  object
}

prepare_existing_reference_object <- function(
  object,
  dataset,
  donor_crosswalk_file = NULL,
  canonical_metadata_file = NULL
) {
  meta <- if (is.null(canonical_metadata_file)) {
    canonicalize_existing_reference_metadata(
      object,
      dataset,
      donor_crosswalk_file
    )
  } else {
    read_existing_canonical_metadata(
      canonical_metadata_file,
      object,
      dataset
    )
  }
  object@meta.data <- meta
  qualify_object_cells(object, dataset)
}

prepare_analysis_reference_object <- function(
  object,
  dataset,
  canonical_metadata_file = NULL,
  reference_component = "reference"
) {
  validate_processed_object(object)
  meta <- if (is.null(canonical_metadata_file)) {
    object[[]]
  } else {
    read_existing_canonical_metadata(
      canonical_metadata_file,
      object,
      dataset
    )
  }
  required <- c(
    "Analysis_Include", "Analysis_Role", "Original_Cell_ID",
    "Original_Donor_ID", "Original_Specimen_ID",
    "Original_Library_ID"
  )
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L) {
    stop(
      dataset,
      " processed object is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  if (any(!meta$Analysis_Include) ||
    any(meta$Analysis_Role != "reference")) {
    stop(dataset, " analysis object contains excluded cells")
  }
  meta$Reference_Component <- reference_component
  object@meta.data <- meta
  qualify_object_cells(object, dataset)
}

align_metadata_columns <- function(metadata_list) {
  all_columns <- unique(unlist(lapply(metadata_list, names)))
  aligned <- lapply(
    metadata_list,
    function(meta) {
      missing <- setdiff(all_columns, names(meta))
      for (column in missing) {
        meta[[column]] <- NA
      }
      meta[, all_columns, drop = FALSE]
    }
  )
  result <- do.call(rbind, aligned)
  rownames(result) <- result$Cells
  result
}

reference_object_audit <- function(objects_list) {
  do.call(
    rbind,
    lapply(
      names(objects_list),
      function(dataset) {
        object <- objects_list[[dataset]]
        meta <- object[[]]
        data.frame(
          Dataset = dataset,
          Reference_Component = paste(
            sort(unique(meta$Reference_Component)),
            collapse = ";"
          ),
          Cells = ncol(object),
          Features = nrow(object),
          Reported_Donors = length(unique(stats::na.omit(
            meta$Global_Donor_ID
          ))),
          Reported_Specimens = length(unique(stats::na.omit(
            meta$Specimen_ID
          ))),
          Reported_Libraries = length(unique(stats::na.omit(
            meta$Library_ID
          ))),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

integration_method_identity <- function() {
  data.frame(
    Output_Name = c(
      "pca_raw", "pca_scvi", "pca_harmony", "pca_rpca",
      "umap_raw", "umap_scvi", "umap_harmony", "umap_rpca"
    ),
    Display_Name = rep(method_levels, 2L),
    Integration_Method = c(
      "Unintegrated PCA", "scVI", "Harmony", "RPCA",
      "Unintegrated PCA UMAP", "scVI UMAP", "Harmony UMAP", "RPCA UMAP"
    ),
    Reduction = c(
      "pca", "integrated.scvi", "integrated.harmony", "integrated.rpca",
      "umap.unintegrated", "umap.scvi", "umap.harmony", "umap.rpca"
    ),
    Space = c(
      "latent", "latent", "latent", "latent",
      "UMAP", "UMAP", "UMAP", "UMAP"
    ),
    stringsAsFactors = FALSE
  )
}

read_integration_reduction_cell_contract <- function(
  file,
  expected_cells = NULL
) {
  if (!file.exists(file)) {
    stop("Integration reduction audit is missing: ", file)
  }
  audit <- read.delim(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required <- c(
    "Output_Name", "Display_Name", "Reduction", "Space", "Cells",
    "Dimensions", "Cell_Order_SHA256", "Finite"
  )
  missing <- setdiff(required, names(audit))
  if (length(missing) > 0L) {
    stop(
      "Integration reduction audit is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  identity <- integration_method_identity()
  if (anyDuplicated(audit$Output_Name) ||
    !setequal(audit$Output_Name, identity$Output_Name)) {
    stop("Integration reduction audit method identities differ")
  }
  audit <- audit[match(identity$Output_Name, audit$Output_Name), , drop = FALSE]
  identity_columns <- c(
    "Output_Name", "Display_Name", "Reduction", "Space"
  )
  if (!identical(
    audit[, identity_columns, drop = FALSE],
    identity[, identity_columns, drop = FALSE]
  )) {
    stop("Integration reduction audit method definitions differ")
  }
  cells <- suppressWarnings(as.numeric(audit$Cells))
  dimensions <- suppressWarnings(as.numeric(audit$Dimensions))
  finite <- tolower(as.character(audit$Finite)) %in% c("true", "t", "1")
  hashes <- tolower(trimws(as.character(audit$Cell_Order_SHA256)))
  if (any(!is.finite(cells)) || any(cells <= 0) ||
    length(unique(cells)) != 1L ||
    any(!is.finite(dimensions)) || any(dimensions <= 0) ||
    any(!finite) ||
    any(!grepl("^[0-9a-f]{64}$", hashes)) ||
    length(unique(hashes)) != 1L) {
    stop("Integration reduction audit failed cell-contract validation")
  }
  if (!is.null(expected_cells) &&
    (!is.numeric(expected_cells) || length(expected_cells) != 1L ||
      is.na(expected_cells) || cells[[1L]] != expected_cells)) {
    stop("Integration reduction audit cell count differs")
  }
  list(
    Cells = cells[[1L]],
    Cell_Order_SHA256 = hashes[[1L]],
    Audit = audit
  )
}

integration_batch_method_contract <- function() {
  methods <- method_levels[-1L]
  batch_models <- c("dataset", "verified_technical")
  method_grid <- expand.grid(
    Batch_Model = batch_models,
    Method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  method_grid <- method_grid[order(
    match(method_grid$Batch_Model, batch_models),
    match(method_grid$Method, methods)
  ), , drop = FALSE]
  batch_interface <- c(
    RPCA = "Seurat layer split / anchor groups",
    Harmony = "Seurat layer split / Harmony group",
    scVI = "SCVI.setup_anndata(batch_key=...)"
  )
  reference <- c(
    RPCA = "https://satijalab.org/seurat/articles/integration_rpca.html",
    Harmony = "https://doi.org/10.1038/s41592-019-0619-0",
    scVI = "https://docs.scvi-tools.org/en/stable/faq.html"
  )
  data.frame(
    Batch_Model = method_grid$Batch_Model,
    Method = method_grid$Method,
    Batch_Interface = unname(batch_interface[method_grid$Method]),
    Required_Metadata_Key = rep("Integration_Batch_ID", nrow(method_grid)),
    Method_Specific_Batch_Key_Allowed = rep(FALSE, nrow(method_grid)),
    Reference = unname(reference[method_grid$Method]),
    stringsAsFactors = FALSE
  )
}

validate_integration_batch_method_contract <- function(contract) {
  required <- c(
    "Batch_Model", "Method", "Batch_Interface", "Required_Metadata_Key",
    "Method_Specific_Batch_Key_Allowed", "Reference"
  )
  missing <- setdiff(required, names(contract))
  if (length(missing) > 0L) {
    stop("integration batch-method contract lacks: ", paste(
      missing,
      collapse = ", "
    ))
  }
  expected_methods <- method_levels[-1L]
  expected_models <- c("dataset", "verified_technical")
  expected_grid <- expand.grid(
    Batch_Model = expected_models,
    Method = expected_methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  expected_grid <- expected_grid[order(
    match(expected_grid$Batch_Model, expected_models),
    match(expected_grid$Method, expected_methods)
  ), , drop = FALSE]
  # Stored evaluation tables may use a different row order; method identity is name-based.
  observed_key <- paste(contract$Batch_Model, contract$Method, sep = ":")
  expected_key <- paste(expected_grid$Batch_Model, expected_grid$Method, sep = ":")
  if (anyDuplicated(observed_key) || !setequal(observed_key, expected_key)) {
    stop("Batch-method evaluation must contain each expected method/model pair once")
  }
  contract <- contract[match(expected_key, observed_key), , drop = FALSE]
  if (!identical(as.character(contract$Method), expected_grid$Method) ||
    !identical(
      as.character(contract$Batch_Model),
      expected_grid$Batch_Model
    ) ||
    any(contract$Required_Metadata_Key != "Integration_Batch_ID") ||
    any(as.logical(contract$Method_Specific_Batch_Key_Allowed))) {
    stop(
      "all correction methods must use each prespecified audited batch model"
    )
  }
  invisible(TRUE)
}

compute_lisi_results <- function(
  embedding_data,
  meta_data,
  method_identity = integration_method_identity(),
  label_column = "Dataset",
  space = "latent",
  perplexity = 30,
  nn_eps = 0,
  lisi_fun = NULL
) {
  required_identity <- c(
    "Output_Name", "Display_Name", "Reduction", "Space"
  )
  missing_identity <- setdiff(required_identity, names(method_identity))
  if (length(missing_identity) > 0L) {
    stop(
      "integration method identity lacks columns: ",
      paste(missing_identity, collapse = ", ")
    )
  }
  selected <- method_identity[
    method_identity$Space == space, ,
    drop = FALSE
  ]
  expected_methods <- method_levels
  selected <- selected[match(expected_methods, selected$Display_Name), , drop = FALSE]
  if (anyNA(selected$Output_Name) ||
    !identical(as.character(selected$Display_Name), expected_methods) ||
    anyDuplicated(selected$Output_Name) ||
    anyDuplicated(selected$Display_Name)) {
    stop(
      "method identity must map Raw, scVI, Harmony and RPCA exactly once in ",
      space,
      " space"
    )
  }
  missing_embeddings <- setdiff(
    as.character(selected$Output_Name),
    names(embedding_data)
  )
  if (length(missing_embeddings) > 0L) {
    stop(
      "integration embeddings are missing: ",
      paste(missing_embeddings, collapse = ", ")
    )
  }
  if (!is.data.frame(meta_data) || !label_column %in% names(meta_data)) {
    stop("LISI metadata lacks label column: ", label_column)
  }
  if (length(perplexity) != 1L ||
    !is.finite(perplexity) ||
    perplexity <= 0) {
    stop("LISI perplexity must be one positive finite value")
  }
  if (length(nn_eps) != 1L ||
    !is.finite(nn_eps) ||
    nn_eps < 0) {
    stop("LISI nearest-neighbor error bound must be non-negative")
  }
  if (is.null(rownames(meta_data)) || anyDuplicated(rownames(meta_data))) {
    stop("LISI metadata requires unique cell row names")
  }
  if (is.null(lisi_fun)) {
    if (!requireNamespace("lisi", quietly = TRUE)) {
      stop("the lisi package is required for integration validation")
    }
    lisi_fun <- lisi::compute_lisi
  }

  results <- vector("list", nrow(selected))
  for (index in seq_len(nrow(selected))) {
    output_name <- selected$Output_Name[[index]]
    display_name <- selected$Display_Name[[index]]
    embedding <- embedding_data[[output_name]]
    if (!is.matrix(embedding) || nrow(embedding) != nrow(meta_data)) {
      stop(output_name, " embedding dimensions differ from LISI metadata")
    }
    if (is.null(rownames(embedding)) ||
      !identical(rownames(embedding), rownames(meta_data))) {
      stop(output_name, " embedding cell order differs from LISI metadata")
    }
    value <- lisi_fun(
      X = embedding,
      meta_data = meta_data,
      label_colnames = label_column,
      perplexity = perplexity,
      nn_eps = nn_eps
    )
    value <- as.data.frame(value)
    if (nrow(value) != nrow(meta_data) || ncol(value) != 1L) {
      stop(output_name, " LISI output must contain one value per cell")
    }
    results[[index]] <- value[[1L]]
    names(results)[[index]] <- display_name
  }
  result <- as.data.frame(results, check.names = FALSE)
  rownames(result) <- rownames(meta_data)
  attr(result, "method_identity") <- selected
  attr(result, "label_column") <- label_column
  attr(result, "space") <- space
  attr(result, "perplexity") <- perplexity
  attr(result, "nn_eps") <- nn_eps
  result
}

lisi_numeric_contract <- function() {
  c(
    source_commit = "a917556310d8d2c66833dcc35aa3d0f4d1b6e0f4",
    tarball_sha256 =
      "b4724c19c1bbcc7891f63ce70cae61da0b8870819981a8a40e864ed1e3caf884",
    build_patch = "CXX11_to_CXX14",
    numeric_patch = "double_precision_probability_normalization_v1",
    numeric_patch_sha256 =
      "6c8d1e25cf62aaa385ef9fff90eeedfb31519768a9853c5d5a116a85837171dd"
  )
}

validate_lisi_numeric_contract <- function(marker_lines) {
  expected <- paste(
    names(lisi_numeric_contract()),
    lisi_numeric_contract(),
    sep = "="
  )
  if (!identical(as.character(marker_lines), unname(expected))) {
    stop(
      "installed LISI does not match the reviewed double-precision contract"
    )
  }
  invisible(TRUE)
}

read_lisi_numeric_contract <- function() {
  if (!requireNamespace("lisi", quietly = TRUE)) {
    stop("the pinned double-precision lisi package is required")
  }
  marker <- file.path(
    find.package("lisi"),
    "BRAINOMICS_PINNED_SOURCE"
  )
  if (!file.exists(marker)) {
    stop("installed LISI lacks its source and numerical provenance marker")
  }
  marker_lines <- readLines(marker, warn = FALSE)
  validate_lisi_numeric_contract(marker_lines)
  lisi_numeric_contract()
}

canonicalize_lisi_lower_bound <- function(values, tolerance = 1e-12) {
  previous_adjustments <- attr(
    values,
    "theoretical_lower_bound_adjustments",
    exact = TRUE
  )
  previous_minimum <- attr(
    values,
    "pre_adjustment_minimum",
    exact = TRUE
  )
  values <- as.numeric(values)
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("LISI lower-bound tolerance must be one non-negative value")
  }
  if (any(!is.finite(values))) {
    stop("LISI values must all be finite")
  }
  minimum_before <- min(values)
  if (minimum_before < 1 - tolerance) {
    stop(
      "LISI value exceeds the permitted numerical lower-bound tolerance: ",
      format(minimum_before, digits = 17)
    )
  }
  if (is.null(previous_adjustments)) {
    previous_adjustments <- 0L
  }
  if (is.null(previous_minimum)) {
    previous_minimum <- minimum_before
  }
  adjusted <- sum(values < 1)
  values[values < 1] <- 1
  attr(values, "theoretical_lower_bound_adjustments") <-
    as.integer(previous_adjustments + adjusted)
  attr(values, "pre_adjustment_minimum") <-
    min(previous_minimum, minimum_before)
  attr(values, "lower_bound_tolerance") <- tolerance
  values
}

write_tsv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(file, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  if (grepl("[.]gz$", file)) {
    connection <- gzfile(temporary, "wt")
    tryCatch(
      write.table(
        x,
        connection,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
      ),
      finally = close(connection)
    )
  } else {
    write.table(
      x,
      temporary,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )
  }
  if (file.exists(file)) {
    unlink(file)
  }
  if (!file.rename(temporary, file)) {
    stop("failed to atomically write table: ", file)
  }
  invisible(file)
}

write_int32_binary <- function(values, file, transform = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  connection <- base::file(file, open = "wb")
  on.exit(close(connection), add = TRUE)
  chunk_size <- 10000000L
  total <- length(values)
  if (total == 0L) {
    return(invisible(file))
  }
  starts <- seq.int(1, total, by = chunk_size)
  for (start in starts) {
    end <- min(total, start + chunk_size - 1L)
    value <- values[start:end]
    if (!is.null(transform)) {
      value <- transform(value)
    }
    value <- as.integer(value)
    writeBin(value, connection, size = 4L, endian = "little")
  }
  invisible(file)
}

sha256_text_vector <- function(values) {
  temporary <- tempfile("brainomics-text-sha256-")
  on.exit(unlink(temporary), add = TRUE)
  writeLines(as.character(values), temporary, useBytes = TRUE)
  processed_file_sha256(temporary)
}

scvi_source_matrix_contract <- function() {
  datasets <- reference_datasets()
  units <- rep("raw non-negative integer gene counts", length(datasets))
  eligible <- rep(TRUE, length(datasets))
  evidence <- c(
    AllenM1 = "Allen Brain Map author count matrix",
    EGAS00001006537 = "author-released five-region raw count matrices",
    GSE104276 = "GEO GSE104276_all_pfc_2394_UMI_count_NOERCC matrix",
    GSE168408 = "author H5AD raw-count layer",
    GSE186538 = "GEO GSE186538_Human_counts.mtx matrix",
    GSE204683 = "author-released raw UMI count matrix",
    GSE207334 = "GEO GSE207334_Multiome_rna_counts.mtx matrix",
    GSE212606 = "GEO EasySci-RNA count matrix",
    GSE217511 = "GEO per-sample 10x count matrices",
    GSE296073 = "author-released RNA count layer",
    GSE67835 = "GEO per-cell HTSeq raw read-count tables",
    GSE81475 = "author-released counts.txt matrix",
    GSE97942 = "author-released regional count matrices",
    HYPOMAP = "author H5AD raw-count layer",
    Li_et_al_2018 = "author fetal count2 and adult umi.raw matrices",
    Ma_et_al_2022 = "author-released RNA count layer",
    PRJCA015229 = "author-released 10x RNA count matrices",
    ROSMAP = "public 10x gene-count matrix",
    SomaMut = "public Seurat RNA counts layer",
    GSE294786 = "GEO GSE294786_all_counts_transposed count matrix",
    Velmeshev_2023 = "CELLxGENE raw-count layer",
    Wang_2025 = "CELLxGENE raw-count layer"
  )
  if (!identical(names(evidence), datasets)) {
    stop("scVI source-matrix contract does not match the reference datasets")
  }
  data.frame(
    Dataset = datasets,
    Source_Matrix_Units = units,
    Source_Provenance_Eligible = eligible,
    Source_Evidence = unname(evidence),
    stringsAsFactors = FALSE
  )
}

audit_sparse_scvi_values <- function(matrix, chunk_size = 10000000L) {
  if (!inherits(matrix, "dgCMatrix")) {
    stop("scVI value audit requires a dgCMatrix")
  }
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("scVI value-audit chunk size must be a positive integer")
  }
  values <- matrix@x
  result <- list(
    Nonzero_Values = as.numeric(length(values)),
    Nonfinite_Values = 0,
    Negative_Values = 0,
    Fractional_Values = 0,
    Values_Over_Int32 = 0,
    Min_Nonzero = Inf,
    Max_Nonzero = -Inf
  )
  if (length(values) == 0L) {
    result$Min_Nonzero <- NA_real_
    result$Max_Nonzero <- NA_real_
    return(result)
  }
  starts <- seq.int(1, length(values), by = chunk_size)
  for (start in starts) {
    end <- min(length(values), start + chunk_size - 1L)
    value <- values[start:end]
    finite <- is.finite(value)
    result$Nonfinite_Values <- result$Nonfinite_Values + sum(!finite)
    value <- value[finite]
    if (length(value) == 0L) {
      next
    }
    result$Negative_Values <- result$Negative_Values + sum(value < 0)
    result$Fractional_Values <-
      result$Fractional_Values + sum(value != floor(value))
    result$Values_Over_Int32 <-
      result$Values_Over_Int32 + sum(value > .Machine$integer.max)
    result$Min_Nonzero <- min(result$Min_Nonzero, min(value))
    result$Max_Nonzero <- max(result$Max_Nonzero, max(value))
  }
  result
}

audit_scvi_input_eligibility <- function(
  objects_list,
  source_contract = scvi_source_matrix_contract()
) {
  if (!is.list(objects_list) ||
    !identical(names(objects_list), reference_datasets())) {
    stop("scVI eligibility audit requires the ordered reference-dataset list")
  }
  required_contract <- c(
    "Dataset", "Source_Matrix_Units", "Source_Provenance_Eligible",
    "Source_Evidence"
  )
  if (!is.data.frame(source_contract) ||
    any(!required_contract %in% names(source_contract)) ||
    !identical(as.character(source_contract$Dataset), names(objects_list)) ||
    anyDuplicated(source_contract$Dataset) ||
    anyNA(source_contract$Source_Provenance_Eligible)) {
    stop("scVI source-matrix contract is incomplete or out of order")
  }
  rows <- vector("list", length(objects_list))
  for (index in seq_along(objects_list)) {
    dataset <- names(objects_list)[[index]]
    object <- objects_list[[index]]
    meta <- object[[]]
    count_layers <- processed_count_layers(object)
    layer_cells <- unlist(lapply(count_layers, colnames), use.names = FALSE)
    if (length(layer_cells) != ncol(object) ||
      anyDuplicated(layer_cells) ||
      !setequal(layer_cells, colnames(object))) {
      stop(dataset, " count layers do not cover every cell exactly once")
    }
    layer_audits <- lapply(count_layers, audit_sparse_scvi_values)
    sum_field <- function(field) {
      sum(vapply(layer_audits, `[[`, numeric(1), field))
    }
    min_values <- vapply(layer_audits, `[[`, numeric(1), "Min_Nonzero")
    max_values <- vapply(layer_audits, `[[`, numeric(1), "Max_Nonzero")
    numeric_eligible <- all(c(
      sum_field("Nonfinite_Values"),
      sum_field("Negative_Values"),
      sum_field("Fractional_Values"),
      sum_field("Values_Over_Int32")
    ) == 0)
    provenance_eligible <- as.logical(
      source_contract$Source_Provenance_Eligible[[index]]
    )
    eligible <- numeric_eligible && provenance_eligible
    reasons <- character()
    if (!numeric_eligible) {
      reasons <- c(reasons, "matrix values are not valid non-negative int32 counts")
    }
    if (!provenance_eligible) {
      reasons <- c(reasons, paste0(
        "source matrix units are ",
        source_contract$Source_Matrix_Units[[index]]
      ))
    }
    metadata_units <- if ("Expression_Units" %in% names(meta)) {
      paste(sort(unique(stats::na.omit(normalize_missing_metadata(
        meta$Expression_Units
      )))), collapse = ";")
    } else {
      NA_character_
    }
    rows[[index]] <- data.frame(
      Dataset = dataset,
      Cells = ncol(object),
      Features = nrow(object),
      Count_Layers = length(count_layers),
      Nonzero_Values = sum_field("Nonzero_Values"),
      Nonfinite_Values = sum_field("Nonfinite_Values"),
      Negative_Values = sum_field("Negative_Values"),
      Fractional_Values = sum_field("Fractional_Values"),
      Values_Over_Int32 = sum_field("Values_Over_Int32"),
      Min_Nonzero = if (all(is.na(min_values))) {
        NA_real_
      } else {
        min(min_values, na.rm = TRUE)
      },
      Max_Nonzero = if (all(is.na(max_values))) {
        NA_real_
      } else {
        max(max_values, na.rm = TRUE)
      },
      Numeric_Count_Eligible = numeric_eligible,
      Source_Matrix_Units = source_contract$Source_Matrix_Units[[index]],
      Metadata_Expression_Units = if (nzchar(metadata_units)) {
        metadata_units
      } else {
        NA_character_
      },
      Source_Provenance_Eligible = provenance_eligible,
      ScVI_Eligible = eligible,
      Exclusion_Reason = if (eligible) {
        NA_character_
      } else {
        paste(reasons, collapse = "; ")
      },
      Source_Evidence = source_contract$Source_Evidence[[index]],
      stringsAsFactors = FALSE
    )
    rm(object, meta, count_layers, layer_audits)
    gc()
  }
  audit <- do.call(rbind, rows)
  if (!identical(as.character(audit$Dataset), reference_datasets()) ||
    sum(audit$Cells) != sum(vapply(objects_list, ncol, numeric(1)))) {
    stop("scVI eligibility audit does not cover the reference list")
  }
  audit
}

export_scvi_input_shards <- function(
  objects_list,
  features,
  output_dir,
  batch_column = "Integration_Batch_ID",
  eligibility_audit,
  expected_cells = NULL,
  overwrite = FALSE
) {
  if (!is.list(objects_list) || is.null(names(objects_list)) ||
    anyDuplicated(names(objects_list)) ||
    any(!names(objects_list) %in% reference_datasets())) {
    stop("scVI export requires an ordered subset of reference datasets")
  }
  required_audit <- c(
    "Dataset", "Source_Matrix_Units", "Source_Provenance_Eligible",
    "Numeric_Count_Eligible", "ScVI_Eligible"
  )
  if (!is.data.frame(eligibility_audit) ||
    any(!required_audit %in% names(eligibility_audit))) {
    stop("scVI export requires a complete eligibility audit")
  }
  eligible_rows <- eligibility_audit[eligibility_audit$ScVI_Eligible, , drop = FALSE]
  if (!identical(as.character(eligible_rows$Dataset), names(objects_list)) ||
    any(!eligible_rows$Source_Provenance_Eligible) ||
    any(!eligible_rows$Numeric_Count_Eligible)) {
    stop("scVI export objects differ from the audited eligible datasets")
  }
  features <- as.character(features)
  if (length(features) < 2L || anyNA(features) || anyDuplicated(features)) {
    stop("scVI export features must be unique and non-missing")
  }
  output_dir <- normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )
  if (dir.exists(output_dir) && !isTRUE(overwrite)) {
    stop("scVI input directory already exists: ", output_dir)
  }
  staging_dir <- paste0(output_dir, ".tmp.", Sys.getpid())
  if (dir.exists(staging_dir)) {
    unlink(staging_dir, recursive = TRUE)
  }
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  installed <- FALSE
  on.exit(
    {
      if (!installed && dir.exists(staging_dir)) {
        unlink(staging_dir, recursive = TRUE)
      }
    },
    add = TRUE
  )
  writeLines(features, file.path(staging_dir, "features.txt"), useBytes = TRUE)

  rows <- vector("list", length(objects_list))
  all_cells <- character()
  for (index in seq_along(objects_list)) {
    dataset <- names(objects_list)[[index]]
    object <- objects_list[[index]]
    object <- SeuratObject::JoinLayers(object)
    meta <- object[[]]
    if (!batch_column %in% names(meta)) {
      stop(dataset, " lacks scVI batch column: ", batch_column)
    }
    batch <- normalize_missing_metadata(meta[[batch_column]])
    if (anyNA(batch)) {
      stop(dataset, " has missing scVI batch values")
    }
    counts <- SeuratObject::LayerData(
      object,
      assay = SeuratObject::DefaultAssay(object),
      layer = "counts"
    )
    if (!inherits(counts, "dgCMatrix")) {
      stop(dataset, " scVI counts are not a dgCMatrix")
    }
    present <- features[features %in% rownames(counts)]
    if (length(present) == 0L) {
      stop(dataset, " has no selected scVI features")
    }
    counts <- counts[present, , drop = FALSE]
    if (!identical(colnames(counts), rownames(meta))) {
      stop(dataset, " scVI matrix and metadata cell order differ")
    }
    if (any(!is.finite(counts@x)) ||
      any(counts@x < 0) ||
      any(counts@x != floor(counts@x)) ||
      any(counts@x > .Machine$integer.max)) {
      stop(dataset, " scVI input is not non-negative integer count data")
    }

    shard_name <- sprintf("%02d_%s", index, dataset)
    shard_dir <- file.path(staging_dir, "shards", shard_name)
    dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)
    data_file <- file.path(shard_dir, "data.int32.bin")
    indices_file <- file.path(shard_dir, "indices.int32.bin")
    indptr_file <- file.path(shard_dir, "indptr.int32.bin")
    cells_file <- file.path(shard_dir, "cells.txt")
    batch_file <- file.path(shard_dir, "batch.txt")

    write_int32_binary(counts@x, data_file)
    feature_index <- match(rownames(counts), features)
    if (anyNA(feature_index)) {
      stop(dataset, " scVI feature index construction failed")
    }
    write_int32_binary(
      counts@i,
      indices_file,
      transform = function(value) feature_index[value + 1L] - 1L
    )
    write_int32_binary(counts@p, indptr_file)
    writeLines(colnames(counts), cells_file, useBytes = TRUE)
    writeLines(batch, batch_file, useBytes = TRUE)

    relative <- function(file) {
      substring(file, nchar(staging_dir) + 2L)
    }
    rows[[index]] <- data.frame(
      Dataset = dataset,
      Order = index,
      Source_Matrix_Units = eligible_rows$Source_Matrix_Units[[index]],
      Source_Provenance_Eligible = TRUE,
      Numeric_Count_Eligible = TRUE,
      ScVI_Eligible = TRUE,
      Cells = ncol(counts),
      Features = length(features),
      Present_Features = nrow(counts),
      Nonzero_Values = length(counts@x),
      Data_File = relative(data_file),
      Data_SHA256 = processed_file_sha256(data_file),
      Indices_File = relative(indices_file),
      Indices_SHA256 = processed_file_sha256(indices_file),
      Indptr_File = relative(indptr_file),
      Indptr_SHA256 = processed_file_sha256(indptr_file),
      Cells_File = relative(cells_file),
      Cells_SHA256 = processed_file_sha256(cells_file),
      Batch_File = relative(batch_file),
      Batch_SHA256 = processed_file_sha256(batch_file),
      stringsAsFactors = FALSE
    )
    all_cells <- c(all_cells, colnames(counts))
    rm(object, counts, meta)
    gc()
  }
  if (anyDuplicated(all_cells)) {
    stop("scVI input contains duplicated cell identifiers")
  }
  if (!is.null(expected_cells) &&
    (length(all_cells) != length(expected_cells) ||
      !setequal(all_cells, expected_cells))) {
    stop("scVI input cell identifiers differ from the merged reference")
  }
  manifest <- do.call(rbind, rows)
  write_tsv(manifest, file.path(staging_dir, "manifest.tsv"))
  audit <- data.frame(
    Datasets = nrow(manifest),
    Reference_Datasets = nrow(eligibility_audit),
    Excluded_Datasets = sum(!eligibility_audit$ScVI_Eligible),
    Cells = sum(manifest$Cells),
    Features = length(features),
    Nonzero_Values = sum(manifest$Nonzero_Values),
    Batch_Column = batch_column,
    Cell_Order_SHA256 = sha256_text_vector(all_cells),
    Feature_Order_SHA256 = sha256_text_vector(features),
    stringsAsFactors = FALSE
  )
  write_tsv(audit, file.path(staging_dir, "input_audit.tsv"))

  backup_dir <- NULL
  if (dir.exists(output_dir)) {
    backup_dir <- paste0(output_dir, ".previous.", Sys.getpid())
    if (!file.rename(output_dir, backup_dir)) {
      stop("failed to stage the previous scVI input directory")
    }
  }
  if (!file.rename(staging_dir, output_dir)) {
    if (!is.null(backup_dir) && dir.exists(backup_dir)) {
      file.rename(backup_dir, output_dir)
    }
    stop("failed to install the scVI input directory")
  }
  installed <- TRUE
  if (!is.null(backup_dir) && dir.exists(backup_dir)) {
    unlink(backup_dir, recursive = TRUE)
  }
  invisible(audit)
}

read_scvi_latent <- function(output_dir, expected_cells, expected_dims = 50L) {
  audit_file <- file.path(output_dir, "scvi_output_audit.tsv")
  cells_file <- file.path(output_dir, "cell_ids.txt")
  latent_file <- file.path(output_dir, "scvi_latent.float32.bin")
  versions_file <- file.path(output_dir, "scvi_versions.tsv")
  required <- c(audit_file, cells_file, latent_file, versions_file)
  missing <- required[!file.exists(required)]
  if (length(missing) > 0L) {
    stop("scVI output is incomplete: ", paste(missing, collapse = ", "))
  }
  audit <- read.delim(
    audit_file,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (nrow(audit) != 1L ||
    as.numeric(audit$Cells[[1L]]) != length(expected_cells) ||
    as.integer(audit$Latent_Dimensions[[1L]]) != expected_dims ||
    as.character(audit$Latent_SHA256[[1L]]) !=
      processed_file_sha256(latent_file)) {
    stop("scVI output audit differs from the expected reference")
  }
  cells <- readLines(cells_file, warn = FALSE)
  if (anyDuplicated(cells) ||
    length(cells) != length(expected_cells) ||
    !setequal(cells, expected_cells)) {
    stop("scVI output cell identifiers differ from the reference")
  }
  connection <- base::file(latent_file, open = "rb")
  on.exit(close(connection), add = TRUE)
  values <- readBin(
    connection,
    what = numeric(),
    n = length(cells) * expected_dims,
    size = 4L,
    endian = "little"
  )
  if (length(values) != length(cells) * expected_dims ||
    any(!is.finite(values))) {
    stop("scVI latent matrix is truncated or non-finite")
  }
  latent <- matrix(
    values,
    nrow = length(cells),
    ncol = expected_dims,
    byrow = TRUE,
    dimnames = list(
      cells,
      paste0("scVI_", seq_len(expected_dims))
    )
  )
  latent[expected_cells, , drop = FALSE]
}
