# Report the project-reader diagnosis selection separately from the later
# processed-input-to-atlas ledger. Stable source evidence is checked before use.
apply_rosmap_source_attrition <- function(dataset, attrition, evidence_dir = "data") {
  filenames <- c(
    crosswalk = "ROSMAP_source_diagnosis_donor_crosswalk.tsv",
    groups = "ROSMAP_source_diagnosis_group_counts.tsv",
    stages = "ROSMAP_source_to_analysis_attrition.tsv"
  )
  evidence <- lapply(filenames, function(name) {
    utils::read.delim(file.path(evidence_dir, name), sep = "\t", quote = "",
      stringsAsFactors = FALSE, check.names = FALSE)
  })
  donors <- evidence$crosswalk
  groups <- evidence$groups
  stages <- evidence$stages
  included <- donors$Diagnosis_RNA_Group_raw == "nonAD" &
    donors$Diagnosis_Pathology_raw == "no"
  nonad <- donors$Diagnosis_RNA_Group_raw == "nonAD"
  valid <- function(condition) {
    if (!isTRUE(condition)) stop("ROSMAP source-selection evidence disagrees with current publication inputs")
  }
  valid(nrow(donors) > 0L && !anyDuplicated(donors$Original_Donor_ID) &&
    all(donors$Dataset == "ROSMAP") && all(is.finite(donors$Source_Cells)) &&
    all(donors$Source_Cells > 0) && all(donors$Source_Cells == as.integer(donors$Source_Cells)) &&
    !anyNA(included) && identical(as.logical(donors$Reader_Include), included))
  valid(all(donors$Reader_Retained_Cells == donors$Source_Cells * included) &&
    all(donors$Reader_Excluded_Cells == donors$Source_Cells * !included))
  source_cells <- sum(donors$Source_Cells)
  retained_cells <- sum(donors$Reader_Retained_Cells)
  excluded_cells <- sum(donors$Reader_Excluded_Cells)
  valid(identical(as.character(stages$Stage), c("1a", "1b", "2")) &&
    all(stages$Dataset == "ROSMAP") && all(groups$Dataset == "ROSMAP"))
  expected_stages <- data.frame(
    Input_Cells = c(source_cells, sum(donors$Source_Cells[nonad]), retained_cells),
    Input_Donors = c(nrow(donors), sum(nonad), sum(included)),
    Retained_Cells = c(sum(donors$Source_Cells[nonad]), retained_cells, retained_cells),
    Retained_Donors = c(sum(nonad), sum(included), sum(included)),
    Excluded_Cells = c(sum(donors$Source_Cells[!nonad]),
      sum(donors$Source_Cells[nonad & !included]), 0),
    Excluded_Donors = c(sum(!nonad), sum(nonad & !included), 0)
  )
  valid(all(as.matrix(stages[names(expected_stages)]) == as.matrix(expected_stages)))
  keys <- paste(donors$Diagnosis_RNA_Group_raw, donors$Diagnosis_Pathology_raw, sep = "::")
  group_keys <- paste(groups$Diagnosis_RNA_Group_raw, groups$Diagnosis_Pathology_raw, sep = "::")
  valid(!anyDuplicated(group_keys) && setequal(keys, group_keys))
  for (i in seq_len(nrow(groups))) {
    selected <- keys == group_keys[[i]]
    valid(groups$Source_Cells[[i]] == sum(donors$Source_Cells[selected]) &&
      groups$Source_Donors[[i]] == sum(selected) &&
      all(as.logical(groups$Reader_Include[[i]]) == included[selected]))
  }
  at <- which(dataset$Dataset == "ROSMAP")
  ar <- which(attrition$Dataset == "ROSMAP")
  valid(length(at) == 1L && length(ar) == 1L)
  valid(dataset$Source_Cells[[at]] == retained_cells &&
    dataset$Analysis_Cells[[at]] == retained_cells &&
    dataset$Known_Donors[[at]] == sum(included) &&
    attrition$Source_Cells[[ar]] == retained_cells &&
    attrition$Analysis_Cells[[ar]] == retained_cells &&
    attrition$Exported_Cells[[ar]] == retained_cells &&
    all(unlist(attrition[ar, grepl("_Excluded_Cells$", names(attrition)), drop = FALSE]) == 0))
  diagnosis <- unique(paste0(
    "ADdiag3types=", donors$Diagnosis_RNA_Group_raw[included],
    "; Pathologic_diagnosis_of_AD=", donors$Diagnosis_Pathology_raw[included]
  ))
  valid(length(diagnosis) == 1L)
  dataset$Diagnosis_Source_Value[at] <- diagnosis
  dataset$Biological_Exclusion_Reason[at] <- paste0(
    "Project reader excluded ", excluded_cells, " of ", source_cells,
    " released RNA nuclei by the joint diagnosis rule before the ", retained_cells,
    "-cell processed integration input; see ROSMAP_source_to_analysis_attrition.tsv ",
    "and ROSMAP_source_diagnosis_donor_crosswalk.tsv"
  )
  attrition$Count_Semantics[ar] <- paste0(
    "Source_Cells is the ", retained_cells, "-cell processed integration input AFTER project-reader ",
    "diagnosis selection, not the complete ", source_cells,
    "-cell released RNA source. This row reports stage 2 only."
  )
  attrition$Attrition_Boundary[ar] <- paste0(
    "No further cells were excluded from the processed integration input to the current atlas. ",
    "The earlier ", excluded_cells, " project-reader diagnosis exclusions are separately recorded in ",
    "ROSMAP_source_to_analysis_attrition.tsv and must not be described as author filtering."
  )
  list(dataset = dataset, attrition = attrition, evidence = evidence, filenames = filenames)
}
