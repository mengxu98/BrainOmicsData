#!/usr/bin/env Rscript
# Export the frozen pipeline results. This script performs no integration or clustering.
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(data.table);library(jsonlite)})
setDTthreads(4L)
a <- commandArgs(TRUE)
stopifnot(length(a)==4L)
frozen <- normalizePath(a[1]); analysis <- normalizePath(a[2]); package <- a[3]; internal <- a[4]
for (d in c(package,internal,file.path(package,c('expression','metadata','embeddings','objects','validation','scripts','provenance')))) dir.create(d,recursive=TRUE,showWarnings=FALSE)
find_unique_file <- function(root, filename) {
 candidates <- list.files(root,recursive=TRUE,full.names=TRUE)
 matches <- candidates[basename(candidates)==filename]
 if(length(matches)!=1L) stop('Expected one ',filename,' under ',root,'; found ',length(matches))
 matches[[1L]]
}
configured_file <- function(variable, root, filename) {
 value <- Sys.getenv(variable,unset='')
 if(nzchar(value)) normalizePath(value,mustWork=TRUE) else find_unique_file(root,filename)
}
log <- function(x) message(format(Sys.time()),' ',x)
json <- function(x,p) write_json(x,p,pretty=TRUE,auto_unbox=TRUE,na='null',digits=16)
known <- function(x) !is.na(x)&nzchar(trimws(as.character(x)))&!tolower(trimws(as.character(x)))%in%c('unknown','not reported','not available','na','nan')
make_ids <- function(x,prefix) {x<-as.character(x);keys<-unique(x[known(x)]);v<-rep(NA_character_,length(x));ok<-known(x);v[ok]<-sprintf(paste0(prefix,'%06d'),match(x[ok],keys));v}
gzwrite <- function(x,path) {
 tmp<-sub('\\.gz$','.partial',path);fwrite(x,tmp,sep='\t',quote=FALSE,na='',nThread=4)
 status<-system2('pigz',c('-n','-f','-p','4',shQuote(tmp)));stopifnot(status==0L)
 stopifnot(file.rename(paste0(tmp,'.gz'),path))
}
# Fail before the expensive frozen-object read when raw-source QC is not ready.
qc_path<-file.path(internal,'percent_mito_by_cell.tsv.gz')
stopifnot(file.exists(qc_path))
qc<-fread(qc_path)
stopifnot(all(c('Cells','Dataset','percent_mito','Source_Status')%in%names(qc)),nrow(qc)==2602031L,!anyDuplicated(qc$Cells),!anyNA(qc$Cells),uniqueN(qc$Dataset)==22L)
stopifnot(all(qc$Source_Status%in%c('COMPUTED_SOURCE_COUNTS','NO_MT_FEATURES','ZERO_TOTAL','SOURCE_CELL_ABSENT')))
stopifnot(all(is.na(qc$percent_mito)|(is.finite(qc$percent_mito)&qc$percent_mito>=0&qc$percent_mito<=100)))
stopifnot(all(!is.na(qc$percent_mito)==(qc$Source_Status=='COMPUTED_SOURCE_COUNTS')))
log('Loading expression object for serialization')
obj <- readRDS(file.path(frozen,'objects_integrated_clustered.rds'))
stopifnot(ncol(obj)==2602031L,nrow(obj)==24659L)
cells<-colnames(obj); genes<-rownames(obj)
stopifnot(!anyDuplicated(cells),!anyDuplicated(genes))
log('Loading age, donor and specimen metadata')
metadata_file<-configured_file('BRAINOMICS_METADATA_FILE',analysis,'metadata_working.rds')
m<-as.data.table(readRDS(metadata_file))
m<-m[match(cells,Cells)]
stopifnot(identical(m$Cells,cells),uniqueN(m$Dataset)==22L)
qc<-qc[match(cells,Cells)]
stopifnot(identical(qc$Cells,cells),identical(qc$Dataset,m$Dataset))
ann<-fread(file.path(internal,'final_cell_annotations.tsv.gz'))
ann<-ann[match(cells,Cells)];stopifnot(identical(ann$Cells,cells),identical(ann$Dataset,m$Dataset),identical(as.character(ann$Cluster),as.character(m$Cluster)),uniqueN(ann$CellType)==12L)
src<-obj[[]]
stopifnot(identical(rownames(src),cells))
raw<-as.character(src$source_cell_type_original_label)
stopifnot(length(raw)==length(cells))
# Preserve the reported original label, including source-specific taxonomy depth.
raw[!known(raw)]<-NA_character_
sidecar_file<-configured_file('BRAINOMICS_SOURCE_ID_SIDECAR',analysis,'cell_source_identifier_sidecar.tsv.gz')
side<-fread(sidecar_file,select=c('Cells','Original_Sample_ID'))
side<-side[match(cells,Cells)];stopifnot(identical(side$Cells,cells))
pub<-sprintf('Cell%07d',seq_along(cells))
sample<-make_ids(m$Specimen_ID_Scoped,'S');donor<-make_ids(m$Canonical_Donor_ID,'D')
libsrc<-as.character(m$Library_ID_Scoped)
libsrc[!grepl('^verified([ :]|$)',tolower(trimws(m$Library_ID_Verification_Status)))]<-NA_character_
lib<-make_ids(libsrc,'L')
ranges<-c('[4,8) PCW','[8,10) PCW','[10,13) PCW','[13,16) PCW','[16,19) PCW','[19,24) PCW','[24,40) PCW','[0,0.5) years','[0.5,1) years','[1,6) years','[6,12) years','[12,20) years','[20,40) years','[40,60) years','[60,+Inf) years')
out<-data.table(Cells=pub,Dataset=m$Dataset,Technology=m$Technology,Sequence=m$Assay_Type,Sample_ID=sample,Original_Sample_ID=side$Original_Sample_ID,Donor_ID=donor,Library_ID=lib,BrainRegion=m$BrainRegion,Age=m$Age,Sex=m$Sex_Source_Standardized,Cluster=as.character(ann$Cluster),seurat_clusters=as.integer(sub('^C','',ann$Cluster)),CellType=ann$CellType,AgeIntervalID=m$AgeIntervalID,AgeInterval=m$AgeInterval,AgeRange=ranges[match(m$AgeIntervalID,paste0('S',1:15))],percent_mito=as.numeric(qc$percent_mito),n_counts=NA_real_,n_genes=NA_integer_,Original_CellType=raw)
stopifnot(ncol(out)==21L,uniqueN(out$Donor_ID)==286L,uniqueN(out$Sample_ID)==448L,uniqueN(out$Library_ID[!is.na(out$Library_ID)])==3744L,!anyNA(out$AgeIntervalID),!anyNA(out$CellType))
gzwrite(data.table(Cells=pub,Internal_Cell_ID=cells,Dataset=m$Dataset,Canonical_Donor_ID=m$Canonical_Donor_ID,Specimen_ID=m$Specimen_ID_Scoped,Library_ID=libsrc),file.path(internal,'private_identifier_crosswalk.tsv.gz'))
origin<-data.table(Dataset=m$Dataset,Original_CellType=raw,Source_Column=as.character(src$source_cell_type_original_label_source_column),Source_Semantics=as.character(src$source_cell_type_original_label_semantics))
# Per-dataset source-column record for Original_CellType; the per-label rows are recoverable from metadata.
semantics_short<-function(x){
  x<-unique(as.character(x))
  out<-vapply(x,function(s){
    if(grepl('retained without harmonization',s,fixed=TRUE)) return('')
    s<-sub('^original ','',s); s<-sub(' annotation$','',s); s<-sub(' final subcluster annotation$',' subcluster',s)
    s
  },character(1))
  paste(sort(unique(out[nzchar(out)])),collapse='; ')
}
ctsrc<-origin[,.(Original_CellType_Source=paste0(paste(sort(unique(Source_Column)),collapse='; '),
  ' (',uniqueN(Original_CellType),' labels)')),by=Dataset]
ctsrc[, Original_CellType_Source := ifelse(nzchar(semantics<-origin[Dataset==.BY[[1L]], semantics_short(Source_Semantics)]),
  paste0(Original_CellType_Source,'; ',semantics), Original_CellType_Source), by=Dataset]
rm(src,ann,side,raw,origin);gc()
layers<-Layers(obj[['RNA']],search='^counts\\.')
stopifnot(length(layers)==22L)
manifest<-list();seen<-rep(FALSE,length(cells))
rm(qc);gc()
for(layer in layers) {
 x<-LayerData(obj,assay='RNA',layer=layer);stopifnot(inherits(x,'dgCMatrix'),identical(rownames(x),genes))
 idx<-match(colnames(x),cells);stopifnot(!anyNA(idx),!any(seen[idx]));seen[idx]<-TRUE
 dataset<-unique(out$Dataset[idx]);stopifnot(length(dataset)==1L)
 log(paste('Exporting counts',dataset,ncol(x),'cells',length(x@x),'nonzero values'))
 dir<-file.path(package,'expression','shards',dataset);dir.create(dir,recursive=TRUE,showWarnings=FALSE)
 out$n_counts[idx]<-Matrix::colSums(x);out$n_genes[idx]<-as.integer(Matrix::colSums(x>0))
 matrixpath<-file.path(dir,'matrix.mtx.gz');mark<-file.path(internal,paste0('counts_',dataset,'.json'))
 if(file.exists(matrixpath)&&file.exists(mark)){
   prior<-read_json(mark,simplifyVector=TRUE)
   stopifnot(prior$cells==ncol(x),prior$genes==nrow(x),prior$nonzero_values==length(x@x))
   log(paste('Reusing completed shard',dataset))
 } else {
   tmp<-file.path(dir,'matrix.mtx.partial')
   writeLines(c('%%MatrixMarket matrix coordinate integer general','% BrainOmicsData retained source counts',paste(nrow(x),ncol(x),length(x@x))),tmp)
   written<-0
   for(lo in seq.int(1L,ncol(x),by=1000L)) {
     hi<-min(lo+999L,ncol(x));start<-x@p[lo]+1;end<-x@p[hi+1L]
     if(end<start)next
     ix<-seq.int(start,end);values<-x@x[ix]
     stopifnot(all(is.finite(values)),all(values>=0),all(values==round(values)))
     dt<-data.table(i=x@i[ix]+1L,j=rep.int(seq.int(lo,hi),diff(x@p[lo:(hi+1L)])),x=values)
     fwrite(dt,tmp,sep=' ',col.names=FALSE,quote=FALSE,append=TRUE,nThread=4)
     written<-written+nrow(dt)
   }
   stopifnot(written==length(x@x))
   stopifnot(system2('pigz',c('-n','-f','-p','4',shQuote(tmp)))==0L,file.rename(paste0(tmp,'.gz'),matrixpath))
   json(list(cells=ncol(x),genes=nrow(x),nonzero_values=written,all_values_finite_nonnegative_integer=TRUE),mark)
 }
 gzwrite(data.table(genes,genes,'Gene Expression'),file.path(dir,'features.tsv.gz'))
 # 10X feature and barcode files have no header.
 for (z in c('features.tsv.gz')) {
   p<-file.path(dir,z);tmp<-sub('\\.gz$','.noheader',p)
   stopifnot(system(paste('gzip -cd',shQuote(p),'| tail -n +2 | pigz -n -p 4 >',shQuote(tmp)))==0L,file.rename(tmp,p))
 }
 tmp<-file.path(dir,'barcodes.partial');writeLines(pub[idx],tmp)
 stopifnot(system2('pigz',c('-n','-f','-p','4',shQuote(tmp)))==0L,file.rename(paste0(tmp,'.gz'),file.path(dir,'barcodes.tsv.gz')))
 manifest[[dataset]]<-data.table(Dataset=dataset,Relative_Directory=paste0('shards/',dataset),Cells=ncol(x),Features=nrow(x),Nonzero_Values=length(x@x))
 rm(x);gc()
}
stopifnot(all(seen),!anyNA(out$n_counts),!anyNA(out$n_genes))
fwrite(rbindlist(manifest),file.path(package,'expression/shard_manifest.tsv'),sep='\t')
gzwrite(out,file.path(package,'metadata/metadata.tsv.gz'))
log('Exporting eight saved reductions without recomputation')
redmap<-c(pca='unintegrated_pca',integrated.scvi='scvi',integrated.harmony='harmony',integrated.rpca='integrated_pca',umap.unintegrated='unintegrated_umap',umap.scvi='scvi_umap',umap.harmony='harmony_umap',umap.rpca='integrated_umap')
reductions<-list()
for (n in names(redmap)) {
 e<-Embeddings(obj[[n]]);stopifnot(identical(rownames(e),cells),ncol(e)==if(startsWith(n,'umap.'))2L else 50L)
 for(j in seq_len(ncol(e)))stopifnot(all(is.finite(e[,j])))
 dt<-as.data.table(e);dt<-cbind(data.table(cell_id=pub),dt)
 gzwrite(dt,file.path(package,'embeddings',paste0(redmap[[n]],'.tsv.gz')))
 rownames(e)<-pub;reductions[[n]]<-CreateDimReducObject(embeddings=e,key=Key(obj[[n]]),assay='RNA')
 log(paste('Saved',redmap[[n]]))
}
# Lightweight plotting object contains no expression matrix and no source-level columns.
dummy<-sparseMatrix(i=integer(),j=integer(),dims=c(1L,length(pub)),dimnames=list('placeholder',pub))
md<-as.data.frame(out);rownames(md)<-pub
plotobj<-CreateSeuratObject(counts=dummy,meta.data=md)
plotobj[['nCount_RNA']]<-NULL;plotobj[['nFeature_RNA']]<-NULL;plotobj[['orig.ident']]<-NULL
for(n in names(reductions))plotobj[[n]]<-reductions[[n]]
stopifnot(identical(names(plotobj[[]]),names(out)))
saveRDS(plotobj,file.path(package,'objects/objects_celltype_plot.rds'),compress=FALSE)
# Core dataset-level table. Reference, licence and mitochondrial-gene columns are merged
# into this file during package assembly.
manifest <- out[,.(Cells=.N,Donors=uniqueN(Donor_ID),Specimens=uniqueN(Sample_ID),Libraries=uniqueN(Library_ID[!is.na(Library_ID)])),by=Dataset]
shared <- out[,.(n=uniqueN(Dataset)),by=Donor_ID][n>1L]
shared[, Datasets := out[Donor_ID==.BY[[1L]], paste(sort(unique(Dataset)), collapse='; ')], by=Donor_ID]
collapse_ids<-function(ids){
  ids<-sort(unique(ids)); if(!length(ids)) return('')
  num<-as.integer(sub('^[A-Za-z]+','',ids))
  groups<-split(ids,cumsum(c(1L,diff(num)!=1L)))
  vapply(groups,function(g) if(length(g)==1L) g[[1L]] else paste0(g[[1L]],'-',g[[length(g)]]),character(1)) |> paste(collapse=', ')
}
manifest[, Cross_Dataset_Donors := vapply(Dataset, function(d) {
  ids <- shared[grepl(d, Datasets, fixed=TRUE), Donor_ID]
  if (!length(ids)) '' else paste(sprintf('%s (%s)', collapse_ids(ids), shared[Donor_ID %in% ids, Datasets]), collapse='; ')
}, '')]
manifest <- merge(manifest, ctsrc, by='Dataset', all.x=TRUE)
fwrite(manifest, file.path(package,'provenance/dataset_manifest.tsv'), sep='\t')
log('Package files written')
