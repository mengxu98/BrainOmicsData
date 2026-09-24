#!/usr/bin/env Rscript
# Map source-object cells onto the released cohort and join the source MT fractions.
#
# Cell identifiers are matched exactly; unresolved or conflicting mappings stop the
# script and leave the lookup unwritten. A documented source-cell absence stays
# blank, never zero.
#
# Usage: percent_mito_map_to_cohort.R ANALYSIS_DIR DATA_ROOT QC_DIR WORK_DIR
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
a<-commandArgs(TRUE);stopifnot(length(a)==4L)
analysis<-a[1];data_root<-a[2];qc<-a[3];work<-a[4]
dir.create(qc,recursive=TRUE,showWarnings=FALSE);dir.create(work,recursive=TRUE,showWarnings=FALSE)
find_unique_file<-function(root,filename){candidates<-list.files(root,recursive=TRUE,full.names=TRUE);matches<-candidates[basename(candidates)==filename];if(length(matches)!=1L)stop('Expected one ',filename,' under ',root,'; found ',length(matches));matches[[1L]]}
metadata_file<-Sys.getenv('BRAINOMICS_METADATA_FILE',unset='')
if(!nzchar(metadata_file))metadata_file<-find_unique_file(analysis,'metadata_working.rds')
m<-as.data.table(readRDS(metadata_file))[,.(Cells,Dataset,Original_Cell_ID)]
stopifnot(nrow(m)==2602031L,!anyDuplicated(m$Cells),uniqueN(m$Dataset)==22L)
cluster_file<-Sys.getenv('BRAINOMICS_CLUSTER_ASSIGNMENTS',unset='')
if(!nzchar(cluster_file)){
 frozen<-Sys.getenv('BRAINOMICS_FROZEN_DIR',unset='')
 if(!nzchar(frozen))stop('Set BRAINOMICS_FROZEN_DIR or BRAINOMICS_CLUSTER_ASSIGNMENTS')
 cluster_file<-file.path(frozen,'cluster_assignments.rds')
}
cl<-readRDS(cluster_file)
m<-m[match(cl$Cell,Cells)];stopifnot(identical(m$Cells,cl$Cell))
feature<-fread(file.path(qc,'feature_inventory.tsv'));result<-list();mapping<-list();checks<-list()
for(ds in unique(m$Dataset)){
 message('Mapping ',ds)
 z<-copy(m[Dataset==ds]);stopifnot(!anyDuplicated(z$Original_Cell_ID))
 p<-file.path(data_root,'processed',ds,'metadata_canonical.tsv.gz')
 canon<-fread(p,select=c('Cells','Original_Cell_ID'))
 stopifnot(!anyDuplicated(canon$Original_Cell_ID))
 canon<-canon[match(z$Original_Cell_ID,Original_Cell_ID)]
 stopifnot(identical(canon$Original_Cell_ID,z$Original_Cell_ID))
 p<-file.path(data_root,'processed',ds,'metadata_raw.tsv.gz')
 fields<-names(fread(p,nrows=0))
 if(all(c('Cells','Source_Object_Cell_ID')%in%fields)){
  raw<-fread(p,select=c('Cells','Source_Object_Cell_ID'));stopifnot(!anyDuplicated(raw$Cells))
  raw<-raw[match(canon$Cells,Cells)]
  stopifnot(identical(raw$Cells,canon$Cells))
 }else if('Original_Cell_ID'%in%fields){
  # H5AD raw metadata is keyed directly by its source obs index.
  raw<-fread(p,select='Original_Cell_ID');stopifnot(!anyDuplicated(raw$Original_Cell_ID))
  raw<-raw[match(z$Original_Cell_ID,Original_Cell_ID)]
  stopifnot(identical(raw$Original_Cell_ID,z$Original_Cell_ID))
  raw[,Source_Object_Cell_ID:=Original_Cell_ID]
 }else stop('No documented exact source ID bridge for ',ds)
 stopifnot(!anyNA(raw$Source_Object_Cell_ID),!anyDuplicated(raw$Source_Object_Cell_ID))
 z[,Source_Object_Cell_ID:=raw$Source_Object_Cell_ID];mapping[[ds]]<-copy(z)
 qp<-file.path(qc,paste0(ds,'_source_qc.tsv.gz'))
 if(ds=='GSE296073')qp<-file.path(qc,'gse296073_10x',paste0(ds,'_source_qc.tsv.gz'))
 if(file.exists(qp)){
  q<-fread(qp);stopifnot(all(q$Dataset==ds),!anyDuplicated(q$Source_Object_Cell_ID))
  idx<-match(z$Source_Object_Cell_ID,q$Source_Object_Cell_ID)
  z[,`:=`(percent_mito=as.numeric(q$percent_mito[idx]),Source_Status=as.character(q$Status[idx]),Source_File=as.character(q$Source_File[idx]),Source_Layer=as.character(q$Source_Layer[idx]))]
  z[is.na(idx),Source_Status:='SOURCE_CELL_UNMATCHED']
 }else{
  f<-feature[Dataset==ds];no_mt<-nrow(f)>0L&&all(f$Status=='NO_MT_FEATURES')&&!ds%in%c('GSE168408','GSE296073')
  z[,`:=`(percent_mito=NA_real_,Source_Status=if(no_mt)'NO_MT_FEATURES'else'PENDING_SOURCE_QC',Source_File=paste(f$Source_File,collapse=';'),Source_Layer=NA_character_)]
 }
 checks[[ds]]<-z[,.(Cells=.N,Nonmissing=sum(!is.na(percent_mito))),by=.(Dataset,Source_Status)]
 result[[ds]]<-z
}
fwrite(rbindlist(mapping),file.path(work,'source_cell_crosswalk.tsv.gz'),sep='\t',na='',quote=FALSE)
v<-rbindlist(result);v<-v[match(m$Cells,Cells)];stopifnot(identical(v$Cells,m$Cells))
coverage<-rbindlist(checks);fwrite(coverage,file.path(work,'percent_mito_mapping_coverage.tsv'),sep='\t',na='',quote=FALSE)
ready<-all(v$Source_Status%in%c('COMPUTED_SOURCE_COUNTS','NO_MT_FEATURES','ZERO_TOTAL','SOURCE_CELL_ABSENT'))
if(!ready){
 print(coverage)
 stop('Source mitochondrial-count inputs are incomplete')
}
stopifnot(all(!is.na(v$percent_mito)==(v$Source_Status=='COMPUTED_SOURCE_COUNTS')))
tmp<-file.path(work,'percent_mito_by_cell.partial.tsv.gz');fwrite(v,tmp,sep='\t',na='',quote=FALSE)
stopifnot(file.rename(tmp,file.path(work,'percent_mito_by_cell.tsv.gz')))
message('Per-cell mitochondrial fractions written for ',nrow(v),' cells')
