normalize_metadata_text <- function(x) {
  if (!is.character(x)) return(x)
  x <- gsub(intToUtf8(0x2013), "-", x, fixed = TRUE)
  x <- gsub(intToUtf8(0x2014), "-", x, fixed = TRUE)
  x <- gsub(intToUtf8(0x2212), "-", x, fixed = TRUE)
  x <- trimws(x)
  x
}

parse_reported_age_fields <- function(age) {
  age <- trimws(as.character(age))
  is_pcw <- grepl("PCW", age, ignore.case = TRUE)
  nums <- regmatches(age, gregexpr("[0-9]+\\.?[0-9]*", age))
  value <- vapply(nums, function(x) {
    if (length(x) == 0) return(NA_real_)
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
  sex <- trimws(as.character(sex))
  ifelse(sex %in% c("Female", "Male"), sex, "Not reported")
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
    "Prefrontal Cortex" = "Prefrontal cortex",
    "Periventricular" = "Periventricular region",
    "Midbrain ventral" = "Ventral midbrain",
    "Whole" = "Whole brain (prenatal)",
    "Cerebellar cortex" = "Cerebellum"
  )
  out <- ifelse(region %in% names(region_map), unname(region_map[region]), region)
  out
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
    "S1" = "4-8 PCW",
    "S2" = "8-10 PCW",
    "S3" = "10-13 PCW",
    "S4" = "13-16 PCW",
    "S5" = "16-19 PCW",
    "S6" = "19-24 PCW",
    "S7" = "24-38 PCW",
    "S8" = "0-0.5 years",
    "S9" = "0.5-1 years",
    "S10" = "1-6 years",
    "S11" = "6-12 years",
    "S12" = "12-20 years",
    "S13" = "20-40 years",
    "S14" = "40-60 years",
    "S15" = "60+ years"
  )

  if ("Stage" %in% names(meta)) {
    meta$AgeIntervalID <- as.character(meta$Stage)
  }
  if ("AgeIntervalID" %in% names(meta)) {
    meta$AgeInterval <- unname(age_interval[as.character(meta$AgeIntervalID)])
    meta$AgeRange <- unname(age_range[as.character(meta$AgeIntervalID)])
    missing_interval <- is.na(meta$AgeInterval) | is.na(meta$AgeRange)
    if (any(missing_interval)) {
      bad <- unique(meta$AgeIntervalID[missing_interval])
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
    age_fields <- parse_reported_age_fields(meta$Age)
    meta$age_raw <- age_fields$age_raw
    meta$age_value <- age_fields$age_value
    meta$age_unit <- age_fields$age_unit
    meta$age_sort <- age_fields$age_sort
  }
  if ("Sex" %in% names(meta)) {
    meta$sex_raw <- as.character(meta$Sex)
    meta$sex_standardized <- standardize_reported_sex(meta$Sex)
  }
  if (!"BrainRegion" %in% names(meta) && "Brain_Region" %in% names(meta)) {
    meta$BrainRegion <- meta$Brain_Region
  }
  if ("BrainRegion" %in% names(meta)) {
    meta$BrainRegion <- standardize_brain_region(meta$BrainRegion)
  }

  meta
}
