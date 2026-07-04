#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

repo_dir <- normalizePath(value_after("--repo-dir", "."), mustWork = TRUE)
object_file <- normalizePath(
  value_after("--object-file", "../../data/BrainOmicsData/integration/objects_celltype_plot.rds"),
  mustWork = TRUE
)
out_file <- value_after("--out-file", "results/sample_schema_audit.tsv")

source(file.path(repo_dir, "functions", "metadata_schema.R"))
source(file.path(repo_dir, "functions", "sample_schema.R"))

message("Loading metadata object: ", object_file)
object <- readRDS(object_file)
if (!inherits(object, "Seurat")) stop("object is not a Seurat object")
meta <- object@meta.data
meta[] <- lapply(meta, as.character)
meta <- add_metadata_schema(meta)
meta <- add_sample_schema(meta)
audit <- sample_schema_audit(meta)
audit <- audit[order(audit$Dataset), ]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.table(audit, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message("Wrote sample schema audit: ", out_file)
message("Global donor count: ", length(unique(meta$Donor_ID)))
message("Global biological sample count: ", length(unique(meta$Sample_ID)))
message("Global library count: ", length(unique(meta$Library_ID)))
