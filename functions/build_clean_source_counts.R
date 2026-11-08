#!/usr/bin/env Rscript
# Build real sparse count chunks while measuring final gene detection. No additional cell QC.
suppressPackageStartupMessages({library(data.table); library(Matrix); library(SeuratObject); library(jsonlite)})
setDTthreads(1)
a <- commandArgs(TRUE)
describe_source <- length(a) == 3L && a[1L] == '--describe-source'
if(describe_source) a <- a[-1L]
stopifnot(length(a) == 2L)
run <- normalizePath(a[1]); ds <- a[2]
input <- file.path(run, 'inputs'); plan <- file.path(input, 'plan')
out <- file.path(run, 'run', if(describe_source) 'source_gene_groups' else 'gene_removal_build', ds)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(!file.exists(file.path(out, 'SUCCESS.json')))
panel <- fread(file.path(plan, 'selected_gene_panel_34387.tsv'))
genes <- panel$Gene_Key; stopifnot(length(genes) == 34387L, !anyDuplicated(genes))
ret <- fread(file.path(plan, 'retained_cells_exploration.tsv.gz'))[Dataset == ds, Cells]
map <- fread(file.path(run, 'code', 'source_row_audit_map.tsv.gz'), na.strings='')[Dataset == ds]
fmt <- fread(file.path(input, 'processed', ds, 'processed_object_format.tsv'))
filename <- if(fmt$Format_Type == 'paired_full_analysis_objects') fmt$Analysis_Object else fmt$Full_Object
message(Sys.time(), ' loading ', ds)
obj <- readRDS(file.path(input, 'processed', ds, filename))
meta <- obj[[]]
original <- if('Original_Cell_ID' %in% names(meta)) as.character(meta$Original_Cell_ID) else colnames(obj)
ids <- ifelse(startsWith(original, paste0(ds, '::')), original, paste0(ds, '::', original))
stopifnot(!anyDuplicated(ids), !anyDuplicated(ret), all(ret %in% ids))
meta <- meta[match(ret, ids), , drop=FALSE]; rownames(meta) <- ret
meta$Dataset <- ds
map[, Object_Feature_ID := Original_Feature_ID]
if(!setequal(map$Object_Feature_ID, rownames(obj))) {
  normalized <- gsub('_', '-', map$Original_Feature_ID, fixed=TRUE)
  stopifnot(!anyDuplicated(normalized), setequal(normalized, rownames(obj)))
  map[, Object_Feature_ID := normalized]
}
assay <- DefaultAssay(obj); layers <- Layers(obj[[assay]], search='^counts')
if(describe_source) {
  native_reports <- group_reports <- list()
  seen <- integer(length(ret))
  for(layer in layers) {
    counts <- LayerData(obj, assay=assay, layer=layer)
    target <- match(ids[match(colnames(counts), colnames(obj))], ret)
    keep <- which(!is.na(target)); if(!length(keep)) next
    lm <- map[match(rownames(counts), Object_Feature_ID)]
    stopifnot(!anyNA(lm$Original_Feature_ID))
    selected <- which(lm$Action == 'candidate')
    indices <- match(lm$Gene_Key[selected], genes)
    stopifnot(!anyNA(indices))
    projection <- sparseMatrix(i=indices, j=seq_along(selected), x=1,
                              dims=c(length(genes), length(selected)))
    detected <- numeric(length(selected)); overlaps <- numeric(length(genes))
    for(start in seq.int(1L, length(keep), by=2000L)) {
      jj <- keep[start:min(length(keep), start+1999L)]
      native <- as(counts[selected, jj, drop=FALSE], 'dgCMatrix')
      stopifnot(all(is.finite(native@x)), all(native@x >= 0), all(native@x == round(native@x)))
      positive <- native > 0
      detected <- detected + Matrix::rowSums(positive)
      overlaps <- overlaps + Matrix::rowSums(projection %*% positive > 1)
      seen[target[jj]] <- seen[target[jj]] + 1L
    }
    current <- data.table(Dataset=ds, Layer=layer,
      Original_Feature_ID=lm$Original_Feature_ID[selected], Gene_Key=lm$Gene_Key[selected],
      Source_Ensembl=lm$Source_Ensembl[selected], Detected_Cells=as.integer(detected))
    groups <- current[, .(Layer_Rows=.N, Positive_Rows=sum(Detected_Cells > 0)), by=Gene_Key][Layer_Rows > 1L]
    groups[, `:=`(Dataset=ds, Layer=layer, Coexpressed_Cells=as.integer(overlaps[match(Gene_Key, genes)]))]
    native_reports[[layer]] <- current
    group_reports[[layer]] <- groups
    rm(native, positive, counts); gc(FALSE)
  }
  stopifnot(all(seen == 1L))
  fwrite(rbindlist(native_reports), file.path(out, 'native_rows.tsv.gz'), sep='\t')
  fwrite(rbindlist(group_reports, fill=TRUE), file.path(out, 'gene_groups.tsv.gz'), sep='\t')
  write_json(list(dataset=ds, retained_cells=length(ret), all_cells_visited_once=TRUE,
                  counts_modified=FALSE),
             file.path(out, 'SUCCESS.json'), pretty=TRUE, auto_unbox=TRUE)
  quit(save='no', status=0L)
}
policy_root <- file.path(run, 'run', 'source_gene_policy', ds)
if(!file.exists(file.path(policy_root, 'applied_gene_row_exclusions.tsv'))) policy_root <- run
evidence <- c(file.path(policy_root, c('applied_gene_row_exclusions.tsv', 'allowed_count_addition.tsv',
                                     'applied_gene_removal_policy.json')),
              normalizePath('functions/build_clean_source_counts.R'))
stopifnot(all(file.exists(evidence)))
hashes <- as.list(tools::md5sum(evidence))
mask <- fread(evidence[1])[Dataset == ds]
allowed <- fread(evidence[2])[Dataset == ds]
reference_root <- Sys.getenv('BRAINOMICS_REFERENCE_RUN_ROOT', unset='')
reference_counts <- NULL
if(nzchar(reference_root)) {
  reference_file <- file.path(reference_root, 'matrices/named24659/sources', ds, 'object.rds')
  stopifnot(file.exists(reference_file))
  reference_object <- readRDS(reference_file)
  reference_counts <- LayerData(reference_object, assay='RNA', layer='counts')
  stopifnot(all(rownames(reference_counts) %in% genes), all(ret %in% colnames(reference_counts)))
}
stopifnot(all(mask$Layer %in% c('*', layers)), all(mask$Original_Feature_ID %in% map$Original_Feature_ID))
seen <- integer(length(ret)); detected <- numeric(length(genes))
before_counts <- after_counts <- numeric(length(ret))
chunks <- coverage <- list(); chunk_index <- 0L
for(layer in layers) {
  counts <- LayerData(obj, assay=assay, layer=layer)
  target <- match(ids[match(colnames(counts), colnames(obj))], ret)
  keep <- which(!is.na(target)); if(!length(keep)) next
  lm <- map[match(rownames(counts), Object_Feature_ID)]
  stopifnot(!anyNA(lm$Original_Feature_ID))
  excluded <- mask[Layer == '*' | Layer == layer, Original_Feature_ID]
  initial <- which(lm$Action == 'candidate')
  selected <- which(lm$Action == 'candidate' & !lm$Original_Feature_ID %in% excluded)
  indices <- match(lm$Gene_Key[selected], genes); stopifnot(!anyNA(indices))
  projection <- sparseMatrix(i=indices, j=seq_along(selected), x=1,
                             dims=c(length(genes), length(selected)), dimnames=list(genes, NULL))
  allow_indices <- match(allowed[Layer == layer, Gene_Key], genes)
  safe_indices <- setdiff(unique(indices), allow_indices)
  representation <- rep('absent_from_source_matrix', length(genes))
  representation[genes %in% lm$Gene_Key] <- 'excluded_all_source_rows'
  representation[unique(indices)] <- 'represented_by_retained_source_rows'
  coverage[[layer]] <- data.table(Dataset=ds, Layer=layer, Gene_Key=genes,
                                  Source_Representation=representation)
  for(start in seq.int(1L, length(keep), by=2000L)) {
    jj <- keep[start:min(length(keep), start+1999L)]; tt <- target[jj]
    native <- as(counts[, jj, drop=FALSE], 'dgCMatrix')
    stopifnot(all(is.finite(native@x)), all(native@x >= 0), all(native@x == round(native@x)))
    selected_counts <- native[selected,,drop=FALSE]
    contributions <- projection %*% (selected_counts > 0)
    # Unsupported positive overlap must never reach output, even if an old audit is stale.
    stopifnot(!any(contributions[safe_indices,,drop=FALSE] > 1))
    cleaned <- as(projection %*% selected_counts, 'dgCMatrix')
    colnames(cleaned) <- ret[tt]
    if(!is.null(reference_counts)) {
      expected <- reference_counts[, match(ret[tt], colnames(reference_counts)), drop=FALSE]
      observed <- cleaned[match(rownames(expected), genes), , drop=FALSE]
      dimnames(observed) <- dimnames(expected)
      difference <- drop0(as(observed - expected, 'dgCMatrix'))
      if(length(difference@x)) stop('BASELINE_COUNT_DIFFERENCE: ', ds, '/', layer,
        '; changed entries=', length(difference@x),
        '; first gene=', rownames(difference)[difference@i[1L] + 1L],
        '; stop before further assembly or integration')
    }
    stopifnot(!anyDuplicated(rownames(cleaned)), all(is.finite(cleaned@x)), all(cleaned@x >= 0))
    before_counts[tt] <- Matrix::colSums(native[initial,,drop=FALSE])
    after_counts[tt] <- Matrix::colSums(cleaned)
    stopifnot(all(after_counts[tt] <= before_counts[tt]),
              isTRUE(all.equal(unname(Matrix::colSums(selected_counts)), unname(after_counts[tt]), tolerance=0)))
    detected <- detected + Matrix::rowSums(cleaned > 0)
    seen[tt] <- seen[tt] + 1L
    chunk_index <- chunk_index + 1L
    chunk_file <- sprintf('counts_%05d.rds', chunk_index)
    saveRDS(cleaned, file.path(out, paste0(chunk_file, '.tmp')), compress=FALSE)
    stopifnot(file.rename(file.path(out, paste0(chunk_file, '.tmp')), file.path(out, chunk_file)))
    chunks[[chunk_index]] <- data.table(File=chunk_file, Layer=layer, Cells=length(tt),
                                       Nonzero=length(cleaned@x), Bytes=file.info(file.path(out, chunk_file))$size,
                                       MD5=unname(tools::md5sum(file.path(out, chunk_file))))
    if(chunk_index %% 10L == 0L) message(Sys.time(), ' wrote ', chunk_index, ' chunks; ', sum(seen), ' cells')
    rm(native, selected_counts, contributions, cleaned)
  }
  rm(counts); gc(FALSE)
}
stopifnot(all(seen == 1L), all(before_counts > 0))
meta$Cleaned_Count_Total <- after_counts
saveRDS(meta, file.path(out, 'cell_metadata.rds'))
fwrite(data.table(Gene_Key=genes, Detected_Cells=detected), file.path(out, 'gene_detection.tsv'), sep='\t')
fwrite(rbindlist(coverage), file.path(out, 'source_feature_representation.tsv.gz'), sep='\t')
fwrite(rbindlist(chunks), file.path(out, 'chunks.tsv'), sep='\t')
fwrite(data.table(Cells=ret, Cleaned_Count_Total=after_counts)[Cleaned_Count_Total == 0], file.path(out, 'zero_cells.tsv'), sep='\t')
stopifnot(identical(hashes, as.list(tools::md5sum(evidence))))
result <- list(dataset=ds, reference_cells=length(ret), cell_membership_modified=FALSE,
               counts_modified=TRUE, all_cells_visited_once=TRUE, chunks=chunk_index,
               candidate_genes=length(genes), detected_genes_10=sum(detected >= 10),
               zero_count_cells=sum(after_counts == 0), before_source_counts=sum(before_counts),
               after_source_counts=sum(after_counts), source_count_loss_fraction=1-sum(after_counts)/sum(before_counts),
               no_unsupported_positive_addition=TRUE,
               retained_final_panel_counts_compared=!is.null(reference_counts),
               output_stage='cleaned sparse counts before cross-source panel assembly')
write_json(result, file.path(out, 'SUCCESS.json'), pretty=TRUE, auto_unbox=TRUE)
message(Sys.time(), ' complete ', ds)
