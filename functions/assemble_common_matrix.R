#!/usr/bin/env Rscript
# Keep all sources as separate sparse count layers; never join a >2^31 nnz matrix.
merge_named_objects <- function(objects_list, genes, cells) {
  stopifnot(!anyDuplicated(cells), length(objects_list) > 1L)
  object <- merge(x=objects_list[[1]], y=objects_list[-1],
                  merge.data=FALSE, merge.dr=FALSE, collapse=FALSE)
  stopifnot(identical(rownames(object), genes), identical(colnames(object), cells))
  object[['RNA']] <- SeuratObject::AddMetaData(object[['RNA']],
    metadata=objects_list[[1]][['RNA']][[]])
  layers <- SeuratObject::Layers(object[['RNA']], search='^counts')
  stopifnot(length(layers) == length(objects_list))
  for (i in seq_along(objects_list)) {
    before <- SeuratObject::LayerData(objects_list[[i]], assay='RNA', layer='counts')
    after <- SeuratObject::LayerData(object, assay='RNA', layer=layers[i])
    stopifnot(identical(dimnames(before), dimnames(after)),
              identical(before@p, after@p), identical(before@i, after@i),
              identical(before@x, after@x))
  }
  object
}

main <- function(run) {
  suppressPackageStartupMessages({library(Seurat); library(jsonlite); library(data.table)})
  Sys.setenv(BRAINOMICS_RUN_ROOT=run)
  source('functions/data_paths.R'); source('functions/dataset_metadata.R')
  source('functions/processed_object.R'); source('functions/integration.R')
  options(Seurat.object.assay.version='v5'); setDTthreads(1)
  marker <- file.path(run, 'matrices/MATRIX_ACCEPTED.json')
  stopifnot(!file.exists(marker))
  verified <- fromJSON(file.path(run, 'matrices/NAMED_TRANSFER_VERIFIED.json'))
  stopifnot(verified$state == 'NAMED_TRANSFER_VERIFIED',
    verified$manifest_sha256 == processed_file_sha256(file.path(run, 'named_transfer_expected.json')))
  policy_file <- file.path(run, 'named_gene_panel_policy.json')
  policy <- fromJSON(policy_file)
  panel <- fread(file.path(run, 'named_gene_panel_24659.tsv'))
  stopifnot(policy$final_genes == 24659)
  datasets <- reference_datasets()
  retained <- fread(file.path(run, 'inputs/plan/retained_cells_exploration.tsv.gz'))
  objects_list <- setNames(lapply(datasets, function(ds) {
    path <- file.path(run, 'matrices/named24659/sources', ds)
    s <- fromJSON(file.path(path, 'SUCCESS.json'))
    stopifnot(s$genes == 24659, s$zero_count_cells == 0,
      file.info(file.path(path,'object.rds'))$size == s$bytes)
    message(Sys.time(), ' loading ', ds)
    object <- load_processed_object(file.path(path, 'object.rds'))
    cells <- retained[Dataset == ds, Cells]
    validate_processed_object(object, expected_features=24659, expected_cells=length(cells))
    stopifnot(identical(rownames(object), panel$Gene_Key),
              identical(colnames(object), cells), all(object$Dataset == ds),
              all(Matrix::colSums(processed_counts(object)) > 0))
    object
  }), datasets)
  cells <- unlist(lapply(objects_list, colnames), use.names=FALSE)
  stopifnot(length(cells) == 2602031)
  out <- brainomics_results_dir(); dir.create(out, recursive=TRUE, showWarnings=FALSE)
  list_file <- file.path(out, 'objects_list_processed.rds')
  # Level-1 gzip fits this independent run within the available HPC storage.
  temp <- paste0(list_file, '.tmp')
  con <- gzfile(temp, 'wb', compression=1)
  tryCatch(serialize(objects_list, con), finally=close(con))
  stopifnot(file.rename(temp, list_file))
  message(Sys.time(), ' merging separate count layers')
  object <- merge_named_objects(objects_list, panel$Gene_Key, cells)
  rm(objects_list, retained); gc(FALSE)
  validate_processed_object(object, expected_features=24659, expected_cells=2602031)
  stopifnot(identical(as.character(object$Integration_Batch_ID), paste0('dataset:', object$Dataset)))
  file <- file.path(out, 'objects_filtered.rds')
  save_processed_object(object, file, compression='gzip', validate_reload=FALSE)
  artifacts <- setNames(lapply(c(list_file, file), function(p)
    list(bytes=unname(file.info(p)$size), sha256=processed_file_sha256(p))),
    basename(c(list_file, file)))
  record <- list(reference_cells=2602031, reference_genes=24659,
    source_identity_resolved=TRUE, count_aggregation_resolved=TRUE,
    zero_count_cells=0, artifacts=artifacts)
  write_json(record, paste0(marker,'.tmp'), auto_unbox=TRUE, pretty=TRUE, digits=NA)
  stopifnot(file.rename(paste0(marker,'.tmp'), marker))
  message(Sys.time(), ' accepted common matrix')
}

if (sys.nframe() == 0L) main(normalizePath(commandArgs(TRUE)[1]))
