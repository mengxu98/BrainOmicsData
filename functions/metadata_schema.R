normalize_metadata_text <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (!is.character(x)) {
    return(x)
  }
  converted <- iconv(x, from = "UTF-8", to = "UTF-8", sub = NA)
  invalid <- !is.na(x) & is.na(converted)
  if (any(invalid)) {
    stop("metadata text contains invalid UTF-8")
  }
  x <- converted
  x <- gsub(intToUtf8(0x00A0), " ", x, fixed = TRUE)
  x <- gsub(intToUtf8(0x2013), "-", x, fixed = TRUE)
  x <- gsub(intToUtf8(0x2014), "-", x, fixed = TRUE)
  x <- gsub(intToUtf8(0x2212), "-", x, fixed = TRUE)
  x <- trimws(x)
  x
}

brainomics_metadata_schema_version <- function() {
  "1.6.1"
}

brainomics_identifier_columns <- function(column_names) {
  explicit <- c(
    "Cells", "Cell", "cellId", "Barcode", "barcodes",
    "Sample", "case", "orig.ident", "Source_Object_Cell_ID"
  )
  id_suffix <- grepl(
    "(^|[._])(?:id|ids)$",
    column_names,
    ignore.case = TRUE,
    perl = TRUE
  )
  column_names[column_names %in% explicit | id_suffix]
}

read_brainomics_table <- function(
  file,
  sep = "\t",
  na.strings = c("", "NA"),
  check.names = FALSE
) {
  if (!file.exists(file)) {
    stop("Table is missing: ", file)
  }
  compressed <- grepl("[.]gz$", file, ignore.case = TRUE)
  header_connection <- if (compressed) {
    gzfile(file, "rt")
  } else {
    base::file(file, "rt")
  }
  header <- tryCatch(
    strsplit(
      readLines(header_connection, n = 1L, warn = FALSE),
      sep,
      fixed = TRUE
    )[[1L]],
    finally = close(header_connection)
  )
  identifier_columns <- brainomics_identifier_columns(header)
  can_use_fread <- requireNamespace("data.table", quietly = TRUE) &&
    (!compressed || requireNamespace("R.utils", quietly = TRUE))
  if (can_use_fread) {
    return(as.data.frame(data.table::fread(
      file,
      sep = sep,
      na.strings = na.strings,
      data.table = FALSE,
      check.names = check.names,
      colClasses = if (length(identifier_columns) > 0L) {
        list(character = identifier_columns)
      } else {
        NULL
      }
    )))
  }

  connection <- if (compressed) {
    gzfile(file, "rt")
  } else {
    base::file(file, "rt")
  }
  column_classes <- rep(NA_character_, length(header))
  names(column_classes) <- header
  column_classes[identifier_columns] <- "character"
  tryCatch(
    read.delim(
      connection,
      sep = sep,
      quote = "",
      na.strings = na.strings,
      stringsAsFactors = FALSE,
      check.names = check.names,
      colClasses = column_classes
    ),
    finally = close(connection)
  )
}

parse_reported_age_fields <- function(age) {
  age <- trimws(as.character(age))
  is_pcw <- grepl("PCW", age, ignore.case = TRUE)
  nums <- regmatches(age, gregexpr("[0-9]+\\.?[0-9]*", age))
  value <- vapply(nums, function(x) {
    if (length(x) == 0) {
      return(NA_real_)
    }
    mean(as.numeric(x), na.rm = TRUE)
  }, numeric(1))
  data.frame(
    age_raw = age,
    age_value = value,
    age_unit = ifelse(is_pcw, "PCW", "years"),
    age_sort = ifelse(is_pcw, value - 1000, value),
    stringsAsFactors = FALSE
  )
}

standardize_sequencing_modality <- function(sequence) {
  sequence <- trimws(as.character(sequence))
  ifelse(
    grepl("^scRNA", sequence, ignore.case = TRUE),
    "scRNA-seq",
    ifelse(grepl("^snRNA", sequence, ignore.case = TRUE), "snRNA-seq", "Other/unknown")
  )
}

standardize_reported_sex <- function(sex) {
  sex <- normalize_metadata_text(as.character(sex))
  normalized <- tolower(sex)
  missing_values <- is.na(normalized) | normalized %in% c(
    "", "na", "n/a", "unknown", "not reported", "not available"
  )
  standardized <- rep("Other/unspecified", length(normalized))
  standardized[missing_values] <- "Not reported"
  standardized[!missing_values & normalized %in% c("f", "female")] <-
    "Female"
  standardized[!missing_values & normalized %in% c("m", "male")] <-
    "Male"
  standardized
}

add_sex_schema <- function(meta) {
  if (!"Sex" %in% names(meta) && !"Sex_Source_Raw" %in% names(meta)) {
    return(meta)
  }
  source_sex <- if ("Sex_Source_Raw" %in% names(meta)) {
    as.character(meta$Sex_Source_Raw)
  } else {
    as.character(meta$Sex)
  }
  source_sex <- normalize_metadata_text(source_sex)
  documented_standardized <- if (
    "Sex_Source_Standardized" %in% names(meta)
  ) {
    standardize_reported_sex(meta$Sex_Source_Standardized)
  } else {
    rep("Not reported", nrow(meta))
  }
  standardized <- standardize_reported_sex(source_sex)
  use_documented_mapping <- documented_standardized %in% c(
    "Female", "Male"
  )
  standardized[use_documented_mapping] <-
    documented_standardized[use_documented_mapping]

  meta$Sex_Source_Raw <- source_sex
  if (!"Sex_Source_Column" %in% names(meta)) {
    meta$Sex_Source_Column <- ifelse(
      is.na(source_sex),
      NA_character_,
      "Sex"
    )
  }
  if (!"Sex_Assignment_Method" %in% names(meta)) {
    meta$Sex_Assignment_Method <- ifelse(
      is.na(source_sex),
      "not reported",
      "source reported"
    )
  }
  meta$sex_raw <- source_sex
  meta$Sex_Source_Standardized <- ifelse(
    standardized %in% c("Female", "Male"),
    standardized,
    NA_character_
  )
  meta$sex_standardized <- standardized
  meta$sex_provenance <- as.character(meta$Sex_Assignment_Method)
  meta$sex_standardization_status <- ifelse(
    standardized == "Not reported",
    "source value missing",
    ifelse(
      standardized == "Other/unspecified",
      "source value retained but not mapped to Female/Male",
      ifelse(
        use_documented_mapping,
        "standardized using documented source coding or sample table",
        "source value standardized without inference"
      )
    )
  )
  meta$Sex <- ifelse(
    standardized %in% c("Female", "Male"),
    standardized,
    NA_character_
  )

  group_conflict <- function(group) {
    group <- normalize_metadata_text(as.character(group))
    result <- rep(NA, length(group))
    known <- !is.na(group)
    if (!any(known)) {
      return(result)
    }
    split_values <- split(meta$Sex[known], group[known])
    conflicts <- vapply(
      split_values,
      function(value) {
        length(unique(stats::na.omit(value))) > 1L
      },
      logical(1)
    )
    result[known] <- unname(conflicts[group[known]])
    result
  }
  meta$Sex_Donor_Conflict <- if ("Global_Donor_ID" %in% names(meta)) {
    group_conflict(meta$Global_Donor_ID)
  } else {
    rep(NA, nrow(meta))
  }
  meta$Sex_Specimen_Conflict <- if ("Specimen_ID" %in% names(meta)) {
    group_conflict(meta$Specimen_ID)
  } else {
    rep(NA, nrow(meta))
  }
  meta
}

standardize_brain_region <- function(region) {
  region <- normalize_metadata_text(as.character(region))
  region_map <- c(
    "Intra-temporal cortex" = "Temporal cortex",
    "cortical plate" = "Cortical plate",
    "Anterior cingulate cortex " = "Anterior cingulate cortex",
    "Ganglionic eminences" = "Ganglionic eminence",
    "ARC" = "Arcuate nucleus",
    "DMH" = "Dorsomedial hypothalamus",
    "Fx/OT/ac" = "Fornix/Optic tract/Anterior commissure",
    "LH" = "Lateral hypothalamus",
    "LPOA" = "Lateral preoptic area",
    "LTN" = "Lateral tuberal nucleus",
    "MAM" = "Mammillary nuclei",
    "ME" = "Median eminence",
    "MPOA" = "Medial preoptic area",
    "Perivent" = "Periventricular region",
    "POA" = "Preoptic area",
    "Thalamaus" = "Thalamus",
    "PVN" = "Paraventricular nucleus",
    "SCN" = "Suprachiasmatic nucleus",
    "SON" = "Supraoptic nucleus",
    "TMN" = "Tuberomammillary nucleus",
    "Vent" = "Ventricle",
    "VMH" = "Ventromedial hypothalamus",
    "BA 9/46" = "Frontal cortex",
    "BM_9/10/46" = "Dorsolateral prefrontal cortex",
    "DFC" = "Dorsolateral prefrontal cortex",
    "NCX" = "Neocortex",
    "Pallium" = "Pallium",
    "pallium" = "Pallium",
    "CBC" = "Cerebellum",
    "FC" = "Frontal cortex",
    "V1C" = "Primary visual cortex",
    "CGE" = "Ganglionic eminence",
    "dEC" = "Entorhinal cortex",
    "EC Stream" = "Entorhinal cortex",
    "LGE" = "Ganglionic eminence",
    "MGE" = "Ganglionic eminence",
    "H29" = "Entorhinal cortex",
    "H31" = "Entorhinal cortex",
    "H33" = "Entorhinal cortex",
    "H37" = "Entorhinal cortex",
    "H39" = "Entorhinal cortex",
    "H46" = "Entorhinal cortex",
    "H48" = "Entorhinal cortex",
    "H71" = "Entorhinal cortex",
    "Cortex" = "Cerebral cortex",
    "Cortex entorhinal" = "Entorhinal cortex",
    "Cortex frontal" = "Frontal cortex",
    "Cortex temporal" = "Temporal cortex",
    "Cortex parietal" = "Parietal cortex",
    "Cortex occipital" = "Occipital cortex",
    "Cortex hemisphere A" = "Cerebral cortex",
    "Cortex hemisphere B" = "Cerebral cortex",
    "Midbrain dorsal" = "Dorsal midbrain",
    "Subcortex" = "Subcortical region",
    "Caudate+Putamen" = "Caudate-putamen",
    "Dorsolateral prefrontal cortex (BA9)" =
      "Dorsolateral prefrontal cortex",
    "Prefrontal Cortex" = "Prefrontal cortex",
    "Periventricular" = "Periventricular region",
    "Midbrain ventral" = "Ventral midbrain",
    "Whole" = "Whole brain (prenatal)",
    "Cerebellar cortex" = "Cerebellum"
  )
  out <- ifelse(region %in% names(region_map), unname(region_map[region]), region)
  # Several public objects use lower-case spellings for the same anatomical
  # concepts used elsewhere in the resource.  Canonicalize case only for this
  # curated set of exact synonyms; retain every source spelling separately in
  # brain_region_raw and brain_region_source_label.
  canonical_case <- c(
    "caudal ganglionic eminence" = "Caudal ganglionic eminence",
    "lateral ganglionic eminence" = "Lateral ganglionic eminence",
    "medial ganglionic eminence" = "Medial ganglionic eminence",
    "cingulate cortex" = "Cingulate cortex",
    "insular cortex" = "Insular cortex",
    "visual cortex" = "Visual cortex",
    "cerebral cortex" = "Cerebral cortex",
    "frontal cortex" = "Frontal cortex",
    "ganglionic eminence" = "Ganglionic eminence",
    "neocortex" = "Neocortex",
    "prefrontal cortex" = "Prefrontal cortex",
    "primary motor cortex" = "Primary motor cortex",
    "forebrain" = "Forebrain",
    "telencephalon" = "Telencephalon",
    "temporal cortex" = "Temporal cortex"
  )
  case_key <- tolower(out)
  canonical_match <- !is.na(case_key) & case_key %in% names(canonical_case)
  out[canonical_match] <- unname(canonical_case[case_key[canonical_match]])
  out
}

curated_brain_region_ontology <- function(dataset, region) {
  mapping <- data.frame(
    Dataset = c(
      rep("*", 45L),
      "HYPOMAP"
    ),
    Region = c(
      "Prefrontal cortex", "Dorsolateral prefrontal cortex",
      "Hippocampus", "Fornix/Optic tract/Anterior commissure",
      "Cerebral cortex", "Cortical plate", "Entorhinal cortex",
      "Primary motor cortex", "Germinal matrix",
      "Periventricular region", "Anterior cingulate cortex",
      "Whole brain (prenatal)", "Frontal cortex",
      "Ganglionic eminence", "Lateral hypothalamus",
      "Primary visual cortex", "Cerebellum", "Mammillary nuclei",
      "Thalamus", "Medial preoptic area",
      "Ventromedial hypothalamus", "Neocortex", "Arcuate nucleus",
      "Preoptic area", "Midbrain", "Lateral preoptic area",
      "Vascular", "Tuberomammillary nucleus", "Medulla",
      "Dorsomedial hypothalamus", "Lateral tuberal nucleus",
      "Subcortical region", "Suprachiasmatic nucleus",
      "Temporal cortex", "Pons", "Hypothalamus", "Ventricle",
      "Paraventricular nucleus", "Diencephalon", "Hindbrain",
      "Basal ganglia", "Supraoptic nucleus", "Median eminence",
      "Ventral midbrain", "Choroid",
      "Periventricular region"
    ),
    Ontology_Label = c(
      "prefrontal cortex", "dorsolateral prefrontal cortex",
      "hippocampal formation", NA, "cerebral cortex",
      "cortical plate", "entorhinal cortex", "primary motor cortex",
      NA, NA, "anterior cingulate cortex", "brain",
      "frontal cortex", "ganglionic eminence",
      "lateral hypothalamic area", "primary visual cortex",
      "cerebellum", "mammillary body",
      "dorsal plus ventral thalamus", "medial preoptic region",
      "ventromedial nucleus of hypothalamus", "neocortex",
      "arcuate nucleus of hypothalamus", "preoptic area", "midbrain",
      "lateral preoptic nucleus", NA, "tuberomammillary nucleus",
      "medulla oblongata", "dorsomedial nucleus of hypothalamus",
      "lateral tuberal nucleus", "cerebral subcortex",
      "suprachiasmatic nucleus", "temporal cortex", "pons",
      "hypothalamus", "brain ventricle",
      "paraventricular nucleus of hypothalamus", "diencephalon",
      "hindbrain", "basal ganglion", "supraoptic nucleus",
      "median eminence of neurohypophysis", "midbrain tegmentum",
      "choroid plexus",
      "periventricular zone of hypothalamus"
    ),
    Ontology_ID = c(
      "UBERON:0000451", "UBERON:0009834", "UBERON:0002421", NA,
      "UBERON:0000956", "UBERON:0005343", "UBERON:0002728",
      "UBERON:0001384", NA, NA, "UBERON:0009835",
      "UBERON:0000955", "UBERON:0001870", "UBERON:0004023",
      "UBERON:0002430", "UBERON:0002436", "UBERON:0002037",
      "UBERON:0002206", "UBERON:0001897", "UBERON:0007769",
      "UBERON:0001935", "UBERON:0001950", "UBERON:0001932",
      "UBERON:0001928", "UBERON:0001891", "UBERON:0001931", NA,
      "UBERON:0001936", "UBERON:0001896", "UBERON:0001934",
      "UBERON:0000435", "UBERON:0000454", "UBERON:0002034",
      "UBERON:0016538", "UBERON:0000988", "UBERON:0001898",
      "UBERON:0004086", "UBERON:0001930", "UBERON:0001894",
      "UBERON:0002028", "UBERON:0002420", "UBERON:0001929",
      "UBERON:0002197", "UBERON:0001943", "UBERON:0001886",
      "UBERON:0002271"
    ),
    Anatomical_Level = c(
      "cortical region", "cortical region", "brain region",
      "composite fiber tracts", "brain region", "developmental layer",
      "cortical region", "cortical region", "developmental zone",
      "region with unresolved anatomical scope", "cortical region",
      "whole organ", "cortical region", "developmental brain region",
      "hypothalamic region", "cortical region", "brain region",
      "hypothalamic nuclear complex", "brain region",
      "hypothalamic region", "hypothalamic nucleus", "brain region",
      "hypothalamic nucleus", "hypothalamic region", "brain division",
      "hypothalamic nucleus", "non-regional source category",
      "hypothalamic nucleus", "brain division", "hypothalamic nucleus",
      "hypothalamic nucleus", "brain region", "hypothalamic nucleus",
      "cortical region", "brain division", "brain region",
      "brain cavity", "hypothalamic nucleus", "brain division",
      "brain division", "brain region", "hypothalamic nucleus",
      "neuroendocrine region", "midbrain region",
      "ventricular tissue", "hypothalamic region"
    ),
    Mapping_Method = c(
      rep("curated exact-label ontology mapping", 2L),
      "curated anatomical synonym mapping",
      "not mapped: composite source region",
      rep("curated exact-label ontology mapping", 4L),
      "not mapped: no exact UBERON equivalent",
      "not mapped: anatomical scope is dataset dependent",
      "curated exact-label ontology mapping",
      "curated whole-organ mapping",
      rep("curated exact-label ontology mapping", 2L),
      "curated anatomical synonym mapping",
      rep("curated exact-label ontology mapping", 2L),
      "curated anatomical synonym mapping",
      "curated anatomical synonym mapping",
      "curated anatomical synonym mapping",
      "curated anatomical synonym mapping",
      "curated exact-label ontology mapping",
      "curated dataset-context ontology mapping",
      "curated exact-label ontology mapping",
      "curated exact-label ontology mapping",
      "curated anatomical synonym mapping",
      "not mapped: non-regional source category",
      "curated exact-label ontology mapping",
      "curated anatomical synonym mapping",
      "curated anatomical synonym mapping",
      "curated exact-label ontology mapping",
      "curated anatomical synonym mapping",
      "curated exact-label ontology mapping",
      "curated exact-label ontology mapping",
      "curated exact-label ontology mapping",
      "curated anatomical synonym mapping",
      "curated dataset-context ontology mapping",
      "curated exact-label ontology mapping",
      "curated exact-label ontology mapping",
      "curated anatomical synonym mapping",
      "curated exact-label ontology mapping",
      "curated anatomical synonym mapping",
      "curated dataset-context ontology mapping",
      "curated anatomical synonym mapping",
      "curated anatomical synonym mapping",
      "curated dataset-context ontology mapping"
    ),
    stringsAsFactors = FALSE
  )
  mapping <- rbind(
    mapping,
    data.frame(
      Dataset = "*",
      Region = "Pallium",
      Ontology_Label = "pallium",
      Ontology_ID = "UBERON:0000203",
      Anatomical_Level = "developmental brain region",
      Mapping_Method = "curated exact-label ontology mapping",
      stringsAsFactors = FALSE
    )
  )
  mapping <- rbind(
    mapping,
    data.frame(
      Dataset = rep("*", 9L),
      Region = c(
        "Forebrain", "Striatum", "Brain", "Dorsal midbrain",
        "Parietal cortex", "Occipital cortex", "Caudate-putamen",
        "Telencephalon", "Head"
      ),
      Ontology_Label = c(
        "forebrain", "striatum", "brain", "midbrain tectum",
        "parietal cortex", "occipital cortex", "caudate-putamen",
        "telencephalon", "head"
      ),
      Ontology_ID = c(
        "UBERON:0001890", "UBERON:0002435", "UBERON:0000955",
        "UBERON:0002314", "UBERON:0016530", "UBERON:0016540",
        "UBERON:0005383", "UBERON:0001893", "UBERON:0000033"
      ),
      Anatomical_Level = c(
        "brain division", "brain region", "whole organ",
        "midbrain region", "cortical region", "cortical region",
        "brain region", "brain division", "body region"
      ),
      Mapping_Method = c(
        rep("curated exact-label ontology mapping", 3L),
        "curated anatomical synonym mapping",
        rep("curated exact-label ontology mapping", 5L)
      ),
      stringsAsFactors = FALSE
    )
  )
  if (nrow(mapping) != 56L) {
    stop("curated brain-region ontology mapping is malformed")
  }
  medium_confidence <- mapping$Region %in% c(
    "Hippocampus", "Mammillary nuclei", "Subcortical region",
    "Ventricle", "Ventral midbrain", "Dorsal midbrain", "Choroid"
  )
  mapping$Mapping_Confidence <- ifelse(
    is.na(mapping$Ontology_ID),
    "not applicable",
    ifelse(medium_confidence, "medium", "high")
  )
  mapping$Mapping_Status <- ifelse(
    is.na(mapping$Ontology_ID),
    "pending ontology mapping",
    "ontology mapped"
  )
  mapping$Mapping_Status[
    mapping$Region == "Fornix/Optic tract/Anterior commissure"
  ] <- "composite region; no single ontology term"
  mapping$Mapping_Status[mapping$Region == "Germinal matrix"] <-
    "pending developmental ontology review"
  mapping$Mapping_Status[
    mapping$Dataset == "*" &
      mapping$Region == "Periventricular region"
  ] <- "pending dataset-specific anatomical review"
  mapping$Mapping_Status[mapping$Region == "Vascular"] <-
    "non-regional source category"

  dataset <- as.character(dataset)
  region <- as.character(region)
  result <- mapping[rep(NA_integer_, length(region)), , drop = FALSE]
  for (scope in c("dataset", "default")) {
    candidates <- if (scope == "dataset") {
      paste(dataset, region, sep = "\r")
    } else {
      paste("*", region, sep = "\r")
    }
    keys <- paste(mapping$Dataset, mapping$Region, sep = "\r")
    matched <- match(candidates, keys)
    fill <- is.na(result$Region) & !is.na(matched)
    result[fill, ] <- mapping[matched[fill], , drop = FALSE]
  }
  rownames(result) <- NULL
  result
}

source_cell_type_column_for_dataset <- function(dataset) {
  source_columns <- c(
    AllenM1 = "Cell.Type",
    EGAD00001006049 = "Cell.Type",
    EGAS00001006537 = "Cell.Type",
    GSE103723 = "cell_type",
    GSE104276 = "cell_types",
    GSE144136 = "Cell.Type",
    GSE168408 = "sub_clust",
    GSE178175 = "Cell.Type",
    GSE186538 = "original_name",
    GSE199762 = "all.exp_type",
    GSE202210 = "Cell.Type",
    GSE204683 = "Cell.type",
    GSE207334 = "subclass",
    GSE212606 = "Cell_type",
    GSE217511 = "celltypes",
    GSE261983 = "not retained; placeholder assigned",
    GSE296073 = "cluster1",
    GSE67835 = "Source_Cell_Type",
    GSE81475 = "CellType",
    GSE97942 = "priCluster",
    HYPOMAP = "celltype_annotation",
    Li_et_al_2018 = "Source_Cell_Type_Fine",
    Ma_et_al_2022 = "subclass",
    Nowakowski_et_al_2017 = "original_name",
    PRJCA015229 = "BigCellType",
    ROSMAP = "Celltype",
    SomaMut = "new_clusters3"
  )
  result <- unname(source_columns[as.character(dataset)])
  result[is.na(result)] <- NA_character_
  result
}

add_source_cell_type_schema <- function(meta) {
  if (!"CellType_raw" %in% names(meta) &&
    !"source_cell_type_label" %in% names(meta)) {
    return(meta)
  }
  source_label <- if ("source_cell_type_label" %in% names(meta)) {
    as.character(meta$source_cell_type_label)
  } else {
    as.character(meta$CellType_raw)
  }
  if (!"source_cell_type_label" %in% names(meta)) {
    meta$source_cell_type_label <- source_label
  }
  if (!"source_cell_type_original_label" %in% names(meta)) {
    meta$source_cell_type_original_label <- source_label
  }
  if (!"source_cell_type_level_1" %in% names(meta)) {
    meta$source_cell_type_level_1 <- NA_character_
  }
  if (!"source_cell_type_level_2" %in% names(meta)) {
    meta$source_cell_type_level_2 <- NA_character_
  }
  if (!"source_cell_type_level_3" %in% names(meta)) {
    meta$source_cell_type_level_3 <- NA_character_
  }
  if (!"source_cluster_id" %in% names(meta)) {
    meta$source_cluster_id <- NA_character_
  }
  exact_source_columns <- c(
    "source_cell_type_original_label_source_column",
    "source_cell_type_label_source_column",
    "source_cell_type_level_1_source_column",
    "source_cell_type_level_2_source_column",
    "source_cell_type_level_3_source_column",
    "source_cluster_id_source_column"
  )
  for (column in exact_source_columns) {
    if (!column %in% names(meta)) {
      meta[[column]] <- NA_character_
    }
  }
  if (all(is.na(meta$source_cell_type_original_label_source_column)) &&
    "source_cell_type_label_source_column" %in% names(meta)) {
    meta$source_cell_type_original_label_source_column <-
      meta$source_cell_type_label_source_column
  }
  if (!"source_cell_type_original_label_semantics" %in% names(meta)) {
    meta$source_cell_type_original_label_semantics <- NA_character_
  }
  missing_original_semantics <- is.na(
    meta$source_cell_type_original_label_semantics
  ) | trimws(meta$source_cell_type_original_label_semantics) == ""
  meta$source_cell_type_original_label_semantics[
    missing_original_semantics &
      !is.na(meta$source_cell_type_original_label)
  ] <- "original-study annotation retained without harmonization"
  if (!"source_cell_type_label_semantics" %in% names(meta)) {
    meta$source_cell_type_label_semantics <- NA_character_
  }
  missing_label_semantics <- is.na(meta$source_cell_type_label_semantics) |
    trimws(meta$source_cell_type_label_semantics) == ""
  meta$source_cell_type_label_semantics[
    missing_label_semantics & !is.na(source_label)
  ] <- "source-study preferred annotation retained without harmonization"
  dataset <- if ("Dataset" %in% names(meta)) {
    as.character(meta$Dataset)
  } else {
    rep(NA_character_, nrow(meta))
  }
  inferred_source_column <- source_cell_type_column_for_dataset(dataset)
  if (!"source_cell_type_source_column" %in% names(meta)) {
    meta$source_cell_type_source_column <- inferred_source_column
  } else {
    source_column <- as.character(meta$source_cell_type_source_column)
    missing_source_column <- is.na(source_column) |
      trimws(source_column) == ""
    source_column[missing_source_column] <-
      inferred_source_column[missing_source_column]
    meta$source_cell_type_source_column <- source_column
  }
  if (!"source_cell_type_source_file" %in% names(meta)) {
    meta$source_cell_type_source_file <-
      "retained processed-object metadata"
  } else {
    source_file <- as.character(meta$source_cell_type_source_file)
    missing_source_file <- is.na(source_file) | trimws(source_file) == ""
    source_file[missing_source_file] <-
      "retained processed-object metadata"
    meta$source_cell_type_source_file <- source_file
  }

  normalized_label <- tolower(trimws(source_label))
  missing_label <- is.na(source_label) | trimws(source_label) == ""
  placeholder_label <- !missing_label & normalized_label %in% c(
    "unknown", "unannoted", "unannotated", "not annotated",
    "not reported", "na", "n/a"
  )
  hierarchy_columns <- intersect(
    c(
      "source_cell_type_original_label",
      "source_cell_type_level_1", "source_cell_type_level_2",
      "source_cell_type_level_3", "source_cluster_id"
    ),
    names(meta)
  )
  hierarchy_available <- rep(FALSE, nrow(meta))
  for (column in hierarchy_columns) {
    value <- as.character(meta[[column]])
    normalized <- tolower(trimws(value))
    hierarchy_available <- hierarchy_available |
      (!is.na(value) & trimws(value) != "" & !normalized %in% c(
        "unknown", "unannoted", "unannotated", "not annotated",
        "not reported", "na", "n/a"
      ))
  }
  annotation_available <- (!missing_label & !placeholder_label) |
    hierarchy_available
  if (!"source_cell_type_annotation_status" %in% names(meta)) {
    meta$source_cell_type_annotation_status <- ifelse(
      !annotation_available,
      "not available",
      ifelse(
        missing_label | placeholder_label,
        "partial source hierarchy retained",
        "reported source label retained"
      )
    )
  }
  if (!"source_cell_type_annotation_scope" %in% names(meta)) {
    meta$source_cell_type_annotation_scope <- ifelse(
      !annotation_available,
      "not applicable",
      paste(
        "cell metadata label retained; original cluster/cell derivation",
        "requires source-file audit"
      )
    )
  }
  if (!"source_cell_type_mapping_method" %in% names(meta)) {
    meta$source_cell_type_mapping_method <- ifelse(
      !annotation_available,
      "no source label available",
      "identity; no harmonized-cell-type mapping applied"
    )
  }
  if (!"source_cell_type_confidence" %in% names(meta)) {
    meta$source_cell_type_confidence <- ifelse(
      !annotation_available,
      NA_character_,
      "not reported in retained field"
    )
  }
  if (!"source_cell_type_ontology_id" %in% names(meta)) {
    meta$source_cell_type_ontology_id <- NA_character_
  }
  if (!"source_cell_type_ontology_label" %in% names(meta)) {
    meta$source_cell_type_ontology_label <- NA_character_
  }
  if (!"source_cell_type_doublet_flag" %in% names(meta)) {
    meta$source_cell_type_doublet_flag <- rep(NA, nrow(meta))
  }
  meta$source_cell_type_available <- annotation_available
  meta
}

add_metadata_schema <- function(meta) {
  meta[] <- lapply(meta, normalize_metadata_text)

  age_interval <- c(
    "S1" = "Embryonic",
    "S2" = "Early fetal",
    "S3" = "Early fetal",
    "S4" = "Early mid-fetal",
    "S5" = "Early mid-fetal",
    "S6" = "Late mid-fetal",
    "S7" = "Late fetal",
    "S8" = "Neonatal and early infancy",
    "S9" = "Late infancy",
    "S10" = "Early childhood",
    "S11" = "Middle and late childhood",
    "S12" = "Adolescence",
    "S13" = "Young adulthood",
    "S14" = "Middle adulthood",
    "S15" = "Late adulthood"
  )
  age_range <- c(
    "S1" = "[4, 8) PCW",
    "S2" = "[8, 10) PCW",
    "S3" = "[10, 13) PCW",
    "S4" = "[13, 16) PCW",
    "S5" = "[16, 19) PCW",
    "S6" = "[19, 24) PCW",
    "S7" = "[24, 40) PCW",
    "S8" = "[0, 0.5) years",
    "S9" = "[0.5, 1) years",
    "S10" = "[1, 6) years",
    "S11" = "[6, 12) years",
    "S12" = "[12, 20) years",
    "S13" = "[20, 40) years",
    "S14" = "[40, 60) years",
    "S15" = "[60, +Inf) years"
  )

  if ("Stage" %in% names(meta)) {
    meta$AgeIntervalID <- as.character(meta$Stage)
  }
  if ("AgeIntervalID" %in% names(meta)) {
    meta$AgeInterval <- unname(age_interval[as.character(meta$AgeIntervalID)])
    meta$AgeRange <- unname(age_range[as.character(meta$AgeIntervalID)])
    unmapped_interval <- !is.na(meta$AgeIntervalID) &
      (is.na(meta$AgeInterval) | is.na(meta$AgeRange))
    if (any(unmapped_interval)) {
      bad <- unique(meta$AgeIntervalID[unmapped_interval])
      stop("unmapped Stage/AgeIntervalID values: ", paste(bad, collapse = ", "))
    }
  }

  if ("Sequence" %in% names(meta)) {
    meta$sequencing_modality_raw <- as.character(meta$Sequence)
    meta$sequencing_modality_standardized <- standardize_sequencing_modality(meta$Sequence)
  }
  if ("Technology" %in% names(meta)) {
    meta$sequencing_technology <- as.character(meta$Technology)
  }
  if ("Age" %in% names(meta)) {
    if (all(c("Age_num", "Unit") %in% names(meta))) {
      meta$age_raw <- as.character(meta$Age)
      meta$age_value <- suppressWarnings(as.numeric(meta$Age_num))
      meta$age_unit <- ifelse(
        as.character(meta$Unit) == "PCW",
        "PCW",
        "years"
      )
      meta$age_unit[is.na(meta$Age_num) | is.na(meta$Unit)] <-
        NA_character_
      meta$age_sort <- ifelse(
        meta$age_unit == "PCW",
        meta$age_value - 1000,
        meta$age_value
      )
    } else {
      age_fields <- parse_reported_age_fields(meta$Age)
      meta$age_raw <- age_fields$age_raw
      meta$age_value <- age_fields$age_value
      meta$age_unit <- age_fields$age_unit
      meta$age_sort <- age_fields$age_sort
    }
  }
  meta <- add_sex_schema(meta)
  if (!"BrainRegion" %in% names(meta) && "Brain_Region" %in% names(meta)) {
    meta$BrainRegion <- meta$Brain_Region
  }
  if ("BrainRegion" %in% names(meta)) {
    region_raw <- normalize_metadata_text(as.character(
      meta$BrainRegion
    ))
    region_standardized <- standardize_brain_region(region_raw)
    meta$brain_region_raw <- region_raw
    if (!"brain_region_source_label" %in% names(meta)) {
      meta$brain_region_source_label <- region_raw
    } else {
      source_label <- as.character(meta$brain_region_source_label)
      missing_source_label <- is.na(source_label) |
        trimws(source_label) == ""
      source_label[missing_source_label] <-
        region_raw[missing_source_label]
      meta$brain_region_source_label <- source_label
    }
    meta$brain_region_standardized <- region_standardized
    meta$brain_region_harmonized <- region_standardized
    if (!"brain_region_ontology_id" %in% names(meta)) {
      meta$brain_region_ontology_id <- NA_character_
    }
    if (!"brain_region_ontology_label" %in% names(meta)) {
      meta$brain_region_ontology_label <- NA_character_
    }
    if (!"brain_region_ontology_source" %in% names(meta)) {
      meta$brain_region_ontology_source <- NA_character_
    }
    if (!"brain_region_anatomical_level" %in% names(meta)) {
      meta$brain_region_anatomical_level <- NA_character_
    }
    dataset <- if ("Dataset" %in% names(meta)) {
      as.character(meta$Dataset)
    } else {
      rep(NA_character_, nrow(meta))
    }
    curated_ontology <- curated_brain_region_ontology(
      dataset,
      region_standardized
    )
    curated_match <- !is.na(curated_ontology$Region)
    curated_mapped <- curated_match &
      !is.na(curated_ontology$Ontology_ID)
    missing_ontology_id <- is.na(meta$brain_region_ontology_id) |
      meta$brain_region_ontology_id == ""
    fill_ontology <- missing_ontology_id & curated_mapped
    meta$brain_region_ontology_id[fill_ontology] <-
      curated_ontology$Ontology_ID[fill_ontology]
    missing_ontology_label <- is.na(meta$brain_region_ontology_label) |
      meta$brain_region_ontology_label == ""
    fill_ontology_label <- missing_ontology_label & curated_mapped
    meta$brain_region_ontology_label[fill_ontology_label] <-
      curated_ontology$Ontology_Label[fill_ontology_label]
    missing_ontology_source <- is.na(meta$brain_region_ontology_source) |
      meta$brain_region_ontology_source == ""
    meta$brain_region_ontology_source[
      missing_ontology_source & curated_mapped
    ] <- "UBERON 2026-04-01"
    missing_anatomical_level <- is.na(meta$brain_region_anatomical_level) |
      meta$brain_region_anatomical_level == ""
    meta$brain_region_anatomical_level[
      missing_anatomical_level & curated_match
    ] <- curated_ontology$Anatomical_Level[
      missing_anatomical_level & curated_match
    ]
    ontology_id <- normalize_metadata_text(as.character(
      meta$brain_region_ontology_id
    ))
    ontology_source <- normalize_metadata_text(as.character(
      meta$brain_region_ontology_source
    ))
    ontology_mapped <- !is.na(ontology_id) & ontology_id != ""
    source_ontology <- ontology_mapped &
      !is.na(ontology_source) &
      grepl("CELLxGENE", ontology_source, ignore.case = TRUE)
    base_mapping_method <- ifelse(
      is.na(region_raw) | region_raw == "",
      "not reported",
      ifelse(
        region_raw == region_standardized,
        "source label retained",
        "curated exact-label crosswalk"
      )
    )
    inferred_mapping_method <- ifelse(
      source_ontology,
      "source ontology mapping retained",
      ifelse(
        curated_match,
        curated_ontology$Mapping_Method,
        ifelse(
          ontology_mapped,
          ifelse(
            region_raw == region_standardized,
            "curated exact-label ontology mapping",
            "curated label and ontology crosswalk"
          ),
          base_mapping_method
        )
      )
    )
    if (!"brain_region_mapping_method" %in% names(meta)) {
      meta$brain_region_mapping_method <- inferred_mapping_method
    } else {
      missing_mapping_method <-
        is.na(meta$brain_region_mapping_method) |
          trimws(as.character(meta$brain_region_mapping_method)) == ""
      meta$brain_region_mapping_method[missing_mapping_method] <-
        inferred_mapping_method[missing_mapping_method]
    }
    inferred_mapping_confidence <- ifelse(
      source_ontology,
      "source reported",
      ifelse(
        curated_match,
        curated_ontology$Mapping_Confidence,
        ifelse(
          ontology_mapped,
          "high",
          ifelse(
            base_mapping_method ==
              "curated exact-label crosswalk",
            "high",
            ifelse(
              base_mapping_method ==
                "source label retained",
              "unreviewed",
              "not applicable"
            )
          )
        )
      )
    )
    if (!"brain_region_mapping_confidence" %in% names(meta)) {
      meta$brain_region_mapping_confidence <- inferred_mapping_confidence
    } else {
      missing_mapping_confidence <-
        is.na(meta$brain_region_mapping_confidence) |
          trimws(as.character(meta$brain_region_mapping_confidence)) == ""
      meta$brain_region_mapping_confidence[missing_mapping_confidence] <-
        inferred_mapping_confidence[missing_mapping_confidence]
    }
    if (!"brain_region_laterality" %in% names(meta)) {
      meta$brain_region_laterality <- NA_character_
    }
    if (!"brain_region_source_column" %in% names(meta)) {
      meta$brain_region_source_column <- "BrainRegion canonical alias"
    }
    meta$brain_region_ontology_mapping_status <- ifelse(
      is.na(region_raw) | region_raw == "",
      "not reported",
      ifelse(
        ontology_mapped,
        "ontology mapped",
        ifelse(
          curated_match,
          curated_ontology$Mapping_Status,
          "pending ontology mapping"
        )
      )
    )
    meta$BrainRegion <- region_standardized
  }

  meta <- add_source_cell_type_schema(meta)
  meta
}
