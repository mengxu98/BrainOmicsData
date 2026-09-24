args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1L) normalizePath(args[[1L]]) else normalizePath(".")
lock_file <- file.path(repo_root, "environment", "r-packages.lock.tsv")
if (!file.exists(lock_file)) {
  stop("Missing R package lock manifest: ", lock_file)
}

lock <- utils::read.delim(lock_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(lock$Package), sum(lock$Package == "R") == 1L)
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
  stop("Unified R environment mismatch:\n", paste(detail, collapse = "\n"))
}

# The source marker distinguishes the pinned numerical implementation from the
# upstream package with the same version number.
marker <- file.path(find.package("lisi"), "BRAINOMICS_PINNED_SOURCE")
expected_marker <- readLines(file.path(repo_root, "environment", "lisi", "LISI_PINNED_SOURCE"))
if (!file.exists(marker) || !identical(readLines(marker), expected_marker)) {
  stop("The installed LISI source/patch marker does not match the repository pin")
}

cat(
  "R 4.5.1 environment is ready (", nrow(packages), " packages)\n",
  sep = ""
)
