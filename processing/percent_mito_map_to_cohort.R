#!/usr/bin/env Rscript
# Map source-object cells onto the released cohort and join the source MT fractions.
#
# Cell identifiers are matched exactly; unresolved or conflicting mappings stop the
# script and leave the lookup unwritten. A documented source-cell absence stays
# blank, never zero.
#
# Usage: percent_mito_map_to_cohort.R ANALYSIS_DIR DATA_ROOT QC_DIR WORK_DIR
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(2L)
a<-commandArgs(TRUE);stopifnot(length(a)==4L)
analysis<-a[1];data_root<-a[2];qc<-a[3];work<-a[4]
dir.create(qc,recursive=TRUE,showWarnings=FALSE);dir.create(work,recursive=TRUE,showWarnings=FALSE)
m<-as.data.table(readRDS(file.path(analysis,'07_downstream/revision_20260918/01_metadata/metadata_working.rds')))[,.(Cells,Dataset,Original_Cell_ID)]
stopifnot(nrow(m)==2602031L,!anyDuplicated(m$Cells),uniqueN(m$Dataset)==22L)
cl<-readRDS(file.path(dirname(analysis),'final_results_20260917/cluster_assignments.rds'))
m<-m[match(cl$Cell,Cells)];stopifnot(identical(m$Cells,cl$Cell))
feature<-fread(file.path(qc,'feature_inventory.tsv'));result<-list();mapping<-list();checks<-list()
for(ds in unique(m$Dataset)){
 message(Sys.time(),' Mapping ',ds)
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
if(ready){
 stopifnot(all(!is.na(v$percent_mito)==(v$Source_Status=='COMPUTED_SOURCE_COUNTS')))
 tmp<-file.path(work,'percent_mito_by_cell.partial.tsv.gz');fwrite(v,tmp,sep='\t',na='',quote=FALSE)
 stopifnot(file.rename(tmp,file.path(work,'percent_mito_by_cell.tsv.gz')))
 }
write_json(list(state=if(ready)'LOOKUP_READY'else'INCOMPLETE',cells=nrow(v),nonmissing=sum(!is.na(v$percent_mito)),coverage=coverage,source_ids_one_to_one=TRUE,zero_imputation=FALSE,computed_from_published_panel=FALSE),file.path(work,'mapping_status.json'),pretty=TRUE,auto_unbox=TRUE)
message(Sys.time(),' MAPPING_',if(ready)'COMPLETE'else'INCOMPLETE');print(coverage)
