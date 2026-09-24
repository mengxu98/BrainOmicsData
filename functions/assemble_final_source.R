#!/usr/bin/env Rscript
# Apply the final panel to already-cleaned counts; preserve the frozen cell set.
suppressPackageStartupMessages({library(data.table); library(Matrix); library(SeuratObject); library(jsonlite)})
setDTthreads(1); options(Seurat.object.assay.version='v5')
a <- commandArgs(TRUE); stopifnot(length(a) == 2L)
run <- normalizePath(a[1]); ds <- a[2]
Sys.setenv(BRAINOMICS_RUN_ROOT=run)
source('functions/data_paths.R'); source('functions/metadata_schema.R')
source('functions/dataset_metadata.R'); source('functions/processed_object.R'); source('functions/integration.R')
src <- file.path(run, 'run/gene_removal_build', ds)
out <- file.path(run, 'matrices/sources', ds)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(!file.exists(file.path(out, 'SUCCESS.json')))
panel_file <- file.path(run, 'run/gene_removal_build/final_gene_panel.tsv')
panel <- fread(panel_file); genes <- panel$Gene_Key
summary <- fromJSON(file.path(run, 'run/gene_removal_build/summary.json'))
stopifnot(length(genes) == summary$final_genes, !anyDuplicated(genes), summary$reference_cells == 2602031)
chunks <- fread(file.path(src, 'chunks.tsv'))
stopifnot(sum(chunks$Nonzero) < .Machine$integer.max)
message(Sys.time(), ' assembling ', ds, ': ', nrow(chunks), ' chunks')
parts <- lapply(chunks$File, function(name) {
  path <- file.path(src, name)
  stopifnot(unname(tools::md5sum(path)) == chunks[File == name, MD5])
  x <- readRDS(path); stopifnot(inherits(x, 'dgCMatrix'), all(genes %in% rownames(x)))
  x[genes,,drop=FALSE]
})
counts <- do.call(cbind, parts)
rm(parts); gc(FALSE)
ret <- fread(file.path(run, 'inputs/plan/retained_cells_exploration.tsv.gz'))[Dataset == ds, Cells]
stopifnot(!anyDuplicated(colnames(counts)), setequal(colnames(counts), ret), identical(rownames(counts), genes))
if(!identical(colnames(counts), ret)) counts <- counts[,ret,drop=FALSE]
totals <- Matrix::colSums(counts)
stopifnot(all(totals > 0), all(is.finite(counts@x)), all(counts@x >= 0), all(counts@x == round(counts@x)))
frozen <- readRDS(file.path(src, 'cell_metadata.rds'))
stopifnot(identical(rownames(frozen), ret), all(totals <= frozen$Cleaned_Count_Total))
cleaned_before_panel_total <- sum(frozen$Cleaned_Count_Total)
fmt <- fread(file.path(run, 'inputs/processed', ds, 'processed_object_format.tsv'))
canonical_file <- file.path(run, 'inputs/processed', ds, fmt$Canonical_Metadata)
meta <- as.data.frame(fread(canonical_file))
original <- as.character(meta$Original_Cell_ID)
qualified <- ifelse(startsWith(original, paste0(ds, '::')), original, paste0(ds, '::', original))
stopifnot(!anyDuplicated(qualified), all(ret %in% qualified))
meta <- meta[match(ret, qualified),,drop=FALSE]; rownames(meta) <- ret
meta$Cells <- meta$Integration_Cell_ID <- ret
meta$Dataset <- ds
if('nCount_RNA' %in% names(meta)) meta$Source_nCount_RNA <- meta$nCount_RNA
if('nFeature_RNA' %in% names(meta)) meta$Source_nFeature_RNA <- meta$nFeature_RNA
meta$nCount_RNA <- totals; meta$nFeature_RNA <- Matrix::colSums(counts > 0)
meta$Reference_Component <- if(ds %in% additional_reference_datasets()) 'additional_reference' else 'existing_reference'
meta <- assign_integration_batch(meta, model='dataset')
stopifnot(identical(rownames(meta), colnames(counts)), all(meta$Dataset == ds))
object <- CreateSeuratObject(counts=counts, assay='RNA', project=ds, meta.data=meta,
                             min.cells=0, min.features=0)
rm(counts, frozen); gc(FALSE)
feature_meta <- as.data.frame(panel[, .(Gene_Key, Symbol, Selected_Gene_Type)])
rownames(feature_meta) <- feature_meta$Gene_Key
object[['RNA']] <- AddMetaData(object[['RNA']], metadata=feature_meta)
object@misc$Reprocessing <- list(
    Gene_Row_Policy='configured_gene_panel',
    Source_Representation_File=file.path(src, 'source_feature_representation.tsv.gz'),
    Zero_Padding_Is_Not_Measured_Expression=TRUE)
validate_processed_object(object, expected_features=length(genes), expected_cells=length(ret))
stopifnot(identical(rownames(object), genes), identical(colnames(object), ret),
          length(Layers(object[['RNA']], search='^counts')) == 1L)
file <- file.path(out, 'object.rds')
save_processed_object(object, file, compression='gzip', validate_reload=identical(ds, 'GSE67835'))
saveRDS(object[[]], file.path(out, 'metadata.rds'))
meta_dt <- as.data.table(object[[]])
for(field in c('Global_Donor_ID','Age','AgeInterval','BrainRegion')) {
  stopifnot(field %in% names(meta_dt))
  tab <- meta_dt[, .(Cells=.N), by=field]; setnames(tab, field, 'Value')
  tab[, Field := field]; fwrite(tab, file.path(out, paste0('retained_by_', field, '.tsv')), sep='\t')
}
write_json(list(dataset=ds, reference_cells=length(ret), genes=length(genes), zero_count_cells=0,
    cell_membership_modified=FALSE,
    bytes=file.info(file)$size,
    count_sum=sum(totals), final_panel_removed_counts=cleaned_before_panel_total-sum(totals)),
    file.path(out,'SUCCESS.json'), pretty=TRUE, auto_unbox=TRUE, digits=NA)
message(Sys.time(), ' complete ', ds)
