args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1L) normalizePath(args[[1L]]) else normalizePath(".")
profile <- if (length(args) >= 2L) args[[2L]] else "integration"
lock_name <- switch(profile,
  integration = "r-packages.lock.tsv",
  preprocessing = "r-preprocessing-packages.lock.tsv",
  `local-validation` = "r-local-validation.lock.tsv",
  stop("Unknown environment profile: ", profile)
)
lock_file <- file.path(repo_root, "environment", lock_name)
if (!file.exists(lock_file)) {
  stop("Missing R package lock manifest: ", lock_file)
}

lock <- utils::read.delim(lock_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(as.character(getRversion()), lock$Version[lock$Package == "R"])) {
  stop(
    "R version mismatch: expected ", lock$Version[lock$Package == "R"],
    "; observed ", as.character(getRversion())
  )
}

packages <- lock[lock$Package != "R", , drop = FALSE]
observed <- vapply(packages$Package, function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion(package))
}, character(1L))
version_equal <- mapply(
  function(current, expected) {
    !is.na(current) && identical(
      base::package_version(current),
      base::package_version(expected)
    )
  },
  observed,
  packages$Version,
  USE.NAMES = FALSE
)
mismatch <- !version_equal
if (any(mismatch)) {
  detail <- paste0(
    packages$Package[mismatch], " expected=", packages$Version[mismatch],
    " observed=", ifelse(is.na(observed[mismatch]), "MISSING", observed[mismatch])
  )
  stop("Frozen R environment mismatch:\n", paste(detail, collapse = "\n"))
}

cat(
  "Frozen ", profile, " R environment verified for ",
  nrow(packages), " packages\n",
  sep = ""
)
