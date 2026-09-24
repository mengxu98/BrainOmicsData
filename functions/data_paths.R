brainomics_repo_root <- function() {
  configured <- Sys.getenv("BRAINOMICS_REPO_ROOT", unset = "")
  if (nzchar(configured)) {
    root <- normalizePath(configured, winslash = "/", mustWork = TRUE)
    if (!file.exists(file.path(root, "functions", "data_paths.R"))) {
      stop("BRAINOMICS_REPO_ROOT does not contain the workflow: ", root)
    }
    return(root)
  }
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
  configured <- Sys.getenv("BRAINOMICS_DATA_ROOT", unset = "")
  if (nzchar(configured)) {
    return(normalizePath(configured, winslash = "/", mustWork = FALSE))
  }
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

brainomics_run_root <- function() {
  normalizePath(
    Sys.getenv("BRAINOMICS_RUN_ROOT", unset = file.path(
      brainomics_repo_root(), "results", "run_root"
    )),
    winslash = "/", mustWork = FALSE
  )
}

brainomics_results_dir <- function() {
  normalizePath(
    Sys.getenv("BRAINOMICS_RESULTS_DIR", unset = file.path(
      brainomics_run_root(), "integration"
    )),
    winslash = "/", mustWork = FALSE
  )
}

brainomics_data_path <- function(...) {
  parts <- list(...)
  # Map the logical integration name to the configured results directory.
  if (length(parts) > 0L && is.character(parts[[1L]]) &&
      length(parts[[1L]]) == 1L && grepl("^integration_25(/|$)", parts[[1L]])) {
    suffix <- sub("^integration_25/?", "", parts[[1L]])
    return(do.call(file.path, c(
      list(brainomics_results_dir()),
      if (nzchar(suffix)) list(suffix) else list(), parts[-1L]
    )))
  }
  do.call(file.path, c(list(brainomics_data_root()), parts))
}
