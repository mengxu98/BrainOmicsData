#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("functions/data_paths.R")
source("functions/processed_object.R")
source("functions/integration.R")

lineage_id <- Sys.getenv("BRAINOMICS_LINEAGE_ID", unset = "")
lineage_dir <- brainomics_lineage_path(lineage_id)
decision_file <- file.path(
  "annotation", paste0(lineage_id, "_lineage_subtype_decisions.tsv")
)
if (!dir.exists(lineage_dir) ||
  !lineage_id %chin% c("excitatory", "inhibitory") ||
  !file.exists(decision_file)) {
  stop(
    "A validated excitatory or inhibitory lineage directory and a complete ",
    "cluster-level broad subclass/state decision table are required"
  )
}
evidence_dir_name <- "subtype_evidence"

paths <- list(
  assignments = file.path(
    lineage_dir, "clustering", "lineage_cluster_assignments.tsv.gz"
  ),
  eligibility = file.path(
    lineage_dir, "markers", "lineage_cluster_subtype_eligibility.tsv"
  ),
  subtype_exclusions = file.path(
    lineage_dir, "markers", "lineage_subtype_inference_exclusions.tsv"
  ),
  evidence = file.path(
    lineage_dir,
    evidence_dir_name,
    "lineage_cluster_evidence_summary.tsv"
  ),
  neighborhoods = file.path(
    lineage_dir,
    "neighborhood_evidence",
    "lineage_cluster_neighborhoods.tsv"
  ),
  neighborhood_exit = file.path(
    lineage_dir, "run", "neighborhood_audit.exit"
  )
)
missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0L ||
  !file.exists(file.path(lineage_dir, "_VALIDATED")) ||
  !identical(
    trimws(readLines(file.path(lineage_dir, "run", "markers.exit"))), "0"
  ) ||
  !identical(
    trimws(readLines(paths$neighborhood_exit)), "0"
  )) {
  stop(
    "Subtype assignment requires validated integration, marker and evidence ",
    "outputs: ", paste(missing, collapse = ", ")
  )
}

output_dir <- file.path(lineage_dir, "subtype_annotation")
if (dir.exists(output_dir) &&
  length(list.files(output_dir, all.files = FALSE)) > 0L) {
  stop("Subtype assignment refuses to overwrite completed outputs")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

assignments <- fread(paths$assignments, sep = "\t")
eligibility <- fread(paths$eligibility, sep = "\t")
subtype_exclusions <- fread(paths$subtype_exclusions, sep = "\t")
evidence <- fread(paths$evidence, sep = "\t")
neighborhoods <- fread(paths$neighborhoods, sep = "\t")
decisions <- fread(
  decision_file,
  sep = "\t", na.strings = c("", "NA"),
  keepLeadingZeros = TRUE
)

required_assignment_columns <- c(
  "Cells", "Global_Cluster", "Main_CellType", "Lineage_Cluster"
)
required_decision_columns <- c(
  "Lineage_Cluster", "Lineage_Subtype", "Assignment_Status",
  "Marker_Evidence", "Context_Evidence", "Stability_Evidence",
  "Exclusion_Evidence", "Nomenclature_DOI"
)
if (any(!required_assignment_columns %in% names(assignments)) ||
  any(!required_decision_columns %in% names(decisions)) ||
  anyDuplicated(assignments$Cells) || anyDuplicated(decisions$Lineage_Cluster) ||
  anyDuplicated(eligibility$Lineage_Cluster) ||
  anyDuplicated(evidence$Lineage_Cluster) ||
  anyDuplicated(neighborhoods[, .(Lineage_Cluster, Rank)])) {
  stop("Lineage assignment or subtype-decision schema changed")
}
if (!all(c(
  "Cluster", "Exclude_From_Subtype_Inference", "Excluded_Cells"
) %in% names(subtype_exclusions)) ||
  anyDuplicated(subtype_exclusions$Cluster)) {
  stop("Subtype-inference exclusion schema changed")
}
excluded_global_clusters <- subtype_exclusions[
  as.logical(Exclude_From_Subtype_Inference), as.character(Cluster)
]

cluster_ids <- sort(unique(assignments$Lineage_Cluster))
if (!identical(sort(decisions$Lineage_Cluster), cluster_ids) ||
  !identical(sort(eligibility$Lineage_Cluster), cluster_ids) ||
  !identical(sort(evidence$Lineage_Cluster), cluster_ids) ||
  !identical(sort(unique(neighborhoods$Lineage_Cluster)), cluster_ids)) {
  stop("Subtype decisions must cover every lineage cluster exactly once")
}
expected_neighbor_ranks <- seq_len(min(5L, length(cluster_ids) - 1L))
neighbor_rank_contract <- neighborhoods[
  , .(Ranks = paste(sort(unique(Rank)), collapse = ",")),
  by = Lineage_Cluster
]
if (nrow(neighborhoods) != length(cluster_ids) * length(expected_neighbor_ranks) ||
  any(
    neighbor_rank_contract$Ranks !=
      paste(expected_neighbor_ranks, collapse = ",")
  )) {
  stop("Neighborhood evidence must provide every expected rank per cluster")
}

allowed_status <- c("supported", "unresolved", "excluded")
if (anyNA(decisions$Assignment_Status) ||
  any(!decisions$Assignment_Status %chin% allowed_status)) {
  stop(
    "Assignment_Status must be supported, unresolved or excluded for every ",
    "lineage cluster"
  )
}

allowed_subtypes <- if (lineage_id == "inhibitory") {
  c(
    "PVALB interneurons", "SST interneurons", "VIP interneurons",
    "SNCG interneurons",
    "LAMP5 interneurons",
    "MEIS2/PBX3 LGE/striatal-like inhibitory neurons",
    "Striatal inhibitory neurons",
    "Immature inhibitory neurons", "Purkinje neurons",
    "Cerebellar inhibitory neurons"
  )
} else {
  c(
    "Upper-layer IT-like neurons",
    "RORB-positive L4/5 IT-like neurons",
    "Layer 5 IT-like neurons",
    "L5 ET-like neurons",
    "Layer 6 IT-like neurons",
    "Layer 6 IT CAR3-like neurons",
    "L6 CT-like neurons", "L5/6 NP-like neurons",
    "Layer 6b neurons", "Hippocampal pyramidal neurons",
    "Dentate granule neurons", "Cerebellar granule neurons",
    "Thalamic excitatory neurons", "Immature excitatory state"
  )
}

supported <- decisions$Assignment_Status == "supported"
unsupported <- !supported
required_evidence_columns <- c(
  "Marker_Evidence", "Context_Evidence", "Stability_Evidence",
  "Exclusion_Evidence", "Nomenclature_DOI"
)
if (anyNA(decisions$Lineage_Subtype[supported]) ||
  any(!decisions$Lineage_Subtype[supported] %chin% allowed_subtypes) ||
  any(!is.na(decisions$Lineage_Subtype[unsupported]))) {
  stop(
    paste(
      "Only supported clusters may receive one of the prespecified broad",
      "subclass/state names"
    )
  )
}
for (column in required_evidence_columns) {
  value <- trimws(as.character(decisions[[column]]))
  if (any(supported & (is.na(value) | !nzchar(value)))) {
    stop("Supported subtype decisions require ", column)
  }
}
if (any(
  supported & !grepl(
    "^10\\.[0-9]{4,9}/", as.character(decisions$Nomenclature_DOI)
  )
)) {
  stop("Supported subtype decisions require a DOI-formatted nomenclature source")
}

decision_eligibility <- merge(
  decisions, eligibility,
  by = "Lineage_Cluster", all.x = TRUE, sort = FALSE
)
resolution_stability_columns <- c(
  "Minimum_Resolution_Primary_Recovery",
  "Mean_Resolution_Primary_Recovery",
  "Minimum_Resolution_Alternate_Purity",
  "Mean_Resolution_Alternate_Purity",
  "Minimum_Resolution_Jaccard",
  "Mean_Resolution_Jaccard"
)
if (any(!resolution_stability_columns %in% names(evidence))) {
  stop("Per-cluster resolution stability evidence is incomplete")
}
decision_eligibility <- merge(
  decision_eligibility,
  evidence[, c("Lineage_Cluster", resolution_stability_columns), with = FALSE],
  by = "Lineage_Cluster", all.x = TRUE, sort = FALSE
)
decision_eligibility <- decision_eligibility[
  match(decisions$Lineage_Cluster, Lineage_Cluster)
]
if (anyNA(decision_eligibility$Eligible_Fraction) ||
  any(
    decision_eligibility$Assignment_Status == "supported" &
      decision_eligibility$Eligible_Fraction < 0.99
  )) {
  stop(
    "Clusters with at least 1% formally excluded subtype-inference cells ",
    "cannot receive a supported subtype name"
  )
}

cluster_cells <- assignments[, .(Cells = .N), by = Lineage_Cluster]
decision_eligibility <- merge(
  decision_eligibility, cluster_cells,
  by = "Lineage_Cluster", all.x = TRUE, sort = FALSE
)
if (anyNA(decision_eligibility$Cells) ||
  !identical(
    as.integer(decision_eligibility$Cells),
    as.integer(decision_eligibility$Full_Cells)
  )) {
  stop("Subtype decision cell counts differ from the complete lineage")
}

cell_annotation <- merge(
  assignments[, ..required_assignment_columns],
  decision_eligibility[, c(
    "Lineage_Cluster", "Lineage_Subtype", "Assignment_Status"
  ), with = FALSE],
  by = "Lineage_Cluster", all.x = TRUE, sort = FALSE
)
cell_annotation[, .Cell_Order := match(Cells, assignments$Cells)]
setorder(cell_annotation, .Cell_Order)
cell_annotation[, .Cell_Order := NULL]
cell_annotation[, `:=`(
  Cluster_Assignment_Status = Assignment_Status,
  Subtype_Inference_Eligible = !Global_Cluster %chin% excluded_global_clusters
)]
cell_annotation[Subtype_Inference_Eligible == FALSE, `:=`(
  Lineage_Subtype = NA_character_,
  Assignment_Status = "excluded"
)]
if (nrow(cell_annotation) != nrow(assignments) ||
  anyDuplicated(cell_annotation$Cells) ||
  !identical(cell_annotation$Cells, assignments$Cells) ||
  anyNA(cell_annotation$Assignment_Status) ||
  any(cell_annotation$Subtype_Inference_Eligible == FALSE &
    !is.na(cell_annotation$Lineage_Subtype))) {
  stop("Cell-level subtype sidecar is incomplete or reordered")
}
observed_excluded_cells <- sum(!cell_annotation$Subtype_Inference_Eligible)
expected_excluded_cells <- sum(subtype_exclusions[
  as.logical(Exclude_From_Subtype_Inference), as.numeric(Excluded_Cells)
])
if (observed_excluded_cells != expected_excluded_cells) {
  stop(
    "Cell-level subtype-inference exclusions differ from the marker contract"
  )
}

counts <- cell_annotation[, .(
  Cells = .N,
  Global_Clusters = uniqueN(Global_Cluster),
  Lineage_Clusters = uniqueN(Lineage_Cluster)
), by = .(Assignment_Status, Lineage_Subtype)]
setorder(counts, Assignment_Status, Lineage_Subtype)

files <- c(
  clusters = file.path(output_dir, "lineage_cluster_subtype_annotation.tsv"),
  cells = file.path(output_dir, "lineage_cell_subtype_annotation.tsv.gz"),
  counts = file.path(output_dir, "lineage_subtype_counts.tsv")
)
fwrite(decision_eligibility, files[["clusters"]], sep = "\t", quote = FALSE)
fwrite(cell_annotation, files[["cells"]], sep = "\t", quote = FALSE)
fwrite(counts, files[["counts"]], sep = "\t", quote = FALSE)

contract_file <- file.path(output_dir, "lineage_subtype_annotation_contract.tsv")
contract <- data.table(
  Field = c(
    "Lineage_ID", "Cells", "Clusters", "Supported_Clusters",
    "Unresolved_Clusters", "Excluded_Clusters", "Assignment_Level",
    "Cell_Level_Program_Assignment", "Global_Clusters_Overwritten",
    "Main_CellType_Overwritten",
    "Subtype_Inference_Excluded_Cells",
    "Excluded_Status_Meaning", "Resolution_Stability_Columns"
  ),
  Value = c(
    lineage_id, nrow(cell_annotation), nrow(decision_eligibility),
    sum(decision_eligibility$Assignment_Status == "supported"),
    sum(decision_eligibility$Assignment_Status == "unresolved"),
    sum(decision_eligibility$Assignment_Status == "excluded"),
    "complete-lineage cluster-level broad subclass/state",
    "FALSE", "FALSE", "FALSE",
    observed_excluded_cells,
    paste(
      "excluded from broad subclass/state naming only; retained in complete",
      "lineage integration and clustering"
    ),
    paste(resolution_stability_columns, collapse = ",")
  )
)
fwrite(contract, contract_file, sep = "\t", quote = FALSE)

manifest_file <- file.path(output_dir, "lineage_subtype_annotation_manifest.tsv")
manifest_inputs <- c(files, contract = contract_file)
manifest <- data.table(
  File = basename(manifest_inputs),
  Size_Bytes = as.character(file.info(manifest_inputs)$size)
)
fwrite(manifest, manifest_file, sep = "\t", quote = FALSE)

message(
  "[lineage-subclass-state] assigned supported cluster-level broad labels for ",
  sum(supported), " of ", nrow(decisions), " ", lineage_id,
  " clusters; complete-cell sidecar retains ", nrow(cell_annotation),
  " cells"
)
