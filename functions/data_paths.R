brainomics_repo_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- unique(c(
    current,
    dirname(current),
    dirname(dirname(current))
  ))
  matches <- candidates[
    file.exists(file.path(candidates, "functions", "data_paths.R")) |
      file.exists(file.path(candidates, "run_pipeline.sh"))
  ]
  if (length(matches) == 0L) {
    stop(
      "Cannot locate the BrainOmicsData repository root. ",
      "Run from the repository root or one of its immediate subdirectories."
    )
  }
  matches[[1L]]
}

brainomics_data_root <- function() {
  repo_root <- brainomics_repo_root()
  data_marker <- "/data/BrainOmicsData/"
  marker_start <- regexpr(data_marker, repo_root, fixed = TRUE)[[1L]]
  if (marker_start > 0L) {
    return(substr(
      repo_root,
      1L,
      marker_start + nchar(data_marker) - 2L
    ))
  }
  normalizePath(
    file.path(repo_root, "../../data/BrainOmicsData"),
    winslash = "/",
    mustWork = FALSE
  )
}

brainomics_data_path <- function(...) {
  file.path(brainomics_data_root(), ...)
}

brainomics_lineage_path <- function(lineage_id, analysis_role = "primary") {
  if (!grepl("^[a-z][a-z0-9_]*$", lineage_id)) {
    stop("Lineage ID must be a lowercase filesystem identifier")
  }
  suffix <- switch(analysis_role,
    primary = "complete_primary",
    complete_sensitivity = "complete_sensitivity",
    excluded_sensitivity = "c65_excluded_sensitivity",
    stop("Unknown lineage analysis role: ", analysis_role)
  )
  brainomics_data_path(
    "integration_25",
    "lineage_analysis_20260825_v2",
    paste(lineage_id, suffix, sep = "_")
  )
}
