source("functions/utils.R")

# A cache can have the right category count and still carry obsolete assignments.
read_celltype_assignments <- function(cells = NULL) {
  data.frame(Cells = c("a", "b"), CellType = c("Excitatory neurons", "CGE-derived inhibitory neurons"))
}
final <- read_celltype_assignments()
stopifnot(!inherits(try(validate_celltype_metadata(final, final$Cells), silent = TRUE), "try-error"))
old <- final
old$CellType <- rev(old$CellType)
stopifnot(inherits(try(validate_celltype_metadata(old, old$Cells), silent = TRUE), "try-error"))
old <- final
old$CellType <- NULL
stopifnot(inherits(try(validate_celltype_metadata(old, old$Cells), silent = TRUE), "try-error"))
for (field in c("Detailed_CellType", "Main_CellType", "CellType_Broad")) {
  old <- final
  old[[field]] <- final$CellType
  stopifnot(inherits(try(validate_celltype_metadata(old, old$Cells), silent = TRUE), "try-error"))
}
cat("PASS: exact CellType required; stale, missing and obsolete fields rejected\n")
