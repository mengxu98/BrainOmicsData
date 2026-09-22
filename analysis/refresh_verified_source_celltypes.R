#!/usr/bin/env Rscript
# Refresh only author cell-type metadata in the two sources with verified joins.
suppressPackageStartupMessages({library(data.table); library(SeuratObject)})
source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/integration.R")
backup <- brainomics_data_path("annotation_history", paste0("source_celltypes_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(backup, recursive = TRUE)
for (dataset in c("GSE204683", "GSE104276")) {
  folder <- brainomics_data_path("processed", dataset)
  format <- fread(file.path(folder, "processed_object_format.tsv"))
  metadata_file <- file.path(folder, format$Canonical_Metadata[[1L]])
  old <- as.data.frame(fread(metadata_file))
  revised <- add_existing_source_cell_type_metadata(old, dataset)
  revised$source_cell_type_available <- !is.na(revised$source_cell_type_label)
  changed <- names(revised)[vapply(names(revised), function(k) !identical(old[[k]], revised[[k]]), logical(1))]
  allowed <- grepl("^source_cell_type_|^source_cluster_id", changed) |
    changed %in% c("Cell.type", "CellType_raw", "Source_CellType", "cell_types",
      "Source_Cell_Type_Broad_Original", "Source_Author_Neuronal_Class")
  stopifnot(all(allowed), identical(old$Cells, revised$Cells), !anyDuplicated(revised$Cells))
  destination <- file.path(backup, dataset)
  dir.create(destination)
  crosswalk <- file.path(folder, "source_cell_type_crosswalk.tsv")
  objects <- unique(na.omit(c(format$Full_Object, format$Analysis_Object)))
  for (file in c(metadata_file, crosswalk, file.path(folder, objects))) {
    if (file.exists(file)) stopifnot(file.copy(file, destination))
  }
  for (name in objects) {
    file <- file.path(folder, name)
    object <- readRDS(file)
    at <- match(colnames(object), revised$Cells)
    stopifnot(!anyNA(at), !anyDuplicated(colnames(object)))
    for (field in changed) object@meta.data[[field]] <- revised[[field]][at]
    # No assay, reduction, cell ID, age or sample field is modified.
    saveRDS(object, paste0(file, ".celltypes-new"))
    stopifnot(file.rename(paste0(file, ".celltypes-new"), file))
    rm(object); gc()
  }
  fwrite(revised, metadata_file, sep = "\t", compress = "gzip", na = "NA")
  fwrite(canonical_metadata_crosswalks(revised)$source_cell_type, crosswalk, sep = "\t", na = "NA")
  message(dataset, ": refreshed ", nrow(revised), " cells; missing labels ",
    sum(is.na(old$source_cell_type_label)), " -> ", sum(is.na(revised$source_cell_type_label)))
}
message("Backup: ", backup)
