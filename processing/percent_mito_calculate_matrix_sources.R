#!/usr/bin/env Rscript
# Mitochondrial fraction for sources supplied as Matrix Market or text matrices.
# All source Gene Expression genes stay in the denominator; only documented
# source-object cells are computed (empty 10x droplets are not exported).
#
# Usage: percent_mito_calculate_matrix_sources.R RAW_ROOT OUTPUT_DIR DATASET
suppressPackageStartupMessages({library(Matrix);library(data.table)})
setDTthreads(2L)
a<-commandArgs(TRUE);stopifnot(length(a)==3L);raw<-a[1];out<-a[2];dataset<-a[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE);scores<-list();inventory<-list()
target<-fread(file.path(dirname(raw),'processed',dataset,'metadata_raw.tsv.gz'),select='Source_Object_Cell_ID')[[1]]
stopifnot(!anyDuplicated(target),!anyNA(target))
calculate<-function(path,genes,cells,keep=seq_along(genes),x=NULL){
 message(Sys.time(),' Reading ',path)
 if(is.null(x)){con<-gzfile(path,'rt');x<-readMM(con);close(con)}
 stopifnot(nrow(x)==length(genes),ncol(x)==length(cells),!anyDuplicated(cells))
 source_cells<-length(cells);selected<-which(cells%in%target)
 if(!length(selected))return(invisible(NULL))
 # Limit to documented source-object cells, not millions of empty 10x droplets.
 # All source genes remain in the denominator.
 x<-x[keep,selected,drop=FALSE];cells<-cells[selected];genes<-genes[keep];mt<-grepl('^MT-',sub('^.*[|]','',genes));stopifnot(any(mt))
 if(inherits(x,'sparseMatrix')){
  values<-x@x
  if(length(values))for(i in seq.int(1L,length(values),by=10000000L)){v<-values[i:min(length(values),i+9999999L)];stopifnot(all(is.finite(v)),all(v>=0),all(v==floor(v)))}
 }else for(i in seq_len(ncol(x))){v<-x[,i];stopifnot(all(is.finite(v)),all(v>=0),all(v==floor(v)))}
 total<-Matrix::colSums(x);numerator<-Matrix::colSums(x[mt,,drop=FALSE]);pct<-ifelse(total>0,100*numerator/total,NA_real_)
 stopifnot(all(is.na(pct)|(is.finite(pct)&pct>=0&pct<=100)))
 scores[[path]]<<-data.table(Dataset=dataset,Source_Object_Cell_ID=cells,percent_mito=pct,MT_Counts=numerator,Total_Counts=total,Source_File=path,Source_Layer='source Gene Expression counts',Status=ifelse(total>0,'COMPUTED_SOURCE_COUNTS','ZERO_TOTAL'))
 inventory[[path]]<<-data.table(Dataset=dataset,Source_File=path,Kind='source_count_matrix',Source_Bytes=file.info(path)$size,Count_Semantics='Unmodified source Gene Expression feature universe; exon/intron scope as supplied, not independently verified',Features=length(genes),MT_Genes=sum(mt),MT_Symbols=paste(sort(unique(genes[mt])),collapse=';'),Status='COMPUTED_SOURCE_COUNTS',Source_Cells=source_cells,Computed_Cells=length(cells))
 message(Sys.time(),' Computed ',length(cells),' source cells; MT=',sum(mt));rm(x);gc()
}
dir<-file.path(raw,dataset)
if(dataset%in%c('GSE217511','PRJCA015229','ROSMAP')){
 fs<-list.files(dir,'features.tsv.gz$',recursive=TRUE,full.names=TRUE)
 if(dataset=='PRJCA015229')fs<-fs[grepl('/HM[^/]+/filtered_feature_bc_matrix/',fs)]
 if(dataset=='ROSMAP')fs<-fs[dirname(fs)==file.path(dir,'RNA')]
 stopifnot(length(fs)>0L)
 for(f in fs){
  feat<-fread(f,header=FALSE);genes<-as.character(feat[[if(ncol(feat)>1L)2L else 1L]])
  keep<-if(ncol(feat)>=3L)which(feat[[3]]=='Gene Expression')else seq_along(genes)
  cells<-fread(sub('features.tsv.gz$','barcodes.tsv.gz',f),header=FALSE)[[1]]
  if(dataset=='GSE217511')cells<-paste0(basename(dirname(f)),'_',cells)
  if(dataset=='PRJCA015229')cells<-paste0(cells,'_',basename(dirname(dirname(f))))
  calculate(sub('features.tsv.gz$','matrix.mtx.gz',f),genes,cells,keep)
 }
}else if(dataset=='GSE186538'){
 genes<-fread(file.path(dir,'rawData/GSE186538_Human_genes.txt.gz'),header=FALSE)[[1]]
 cells<-fread(file.path(dir,'rawData/GSE186538_Human_cell_meta.txt.gz'))[[1]]
 calculate(file.path(dir,'rawData/GSE186538_Human_counts.mtx.gz'),genes,cells)
}else if(dataset=='GSE207334'){
 genes<-fread(file.path(dir,'GSE207334_Multiome_rna_genes.txt.gz'),header=FALSE)[[1]]
 cells<-fread(file.path(dir,'GSE207334_Multiome_cell_meta.txt.gz'))[[1]]
 calculate(file.path(dir,'GSE207334_Multiome_rna_counts.mtx.gz'),genes,cells)
}else if(dataset=='GSE212606'){
 genes<-fread(file.path(dir,'GSM6657986_gene_annotation.csv'))$gene_short_name
 cells<-fread(file.path(dir,'GSM6657986_cell_annotation.csv'))[[1]]
 calculate(file.path(dir,'GSM6657986_gene_count.txt.gz'),genes,cells)
}else if(dataset=='GSE81475'){
 path<-file.path(dir,'counts.txt');dt<-fread(path);genes<-dt[[1]];cells<-names(dt)[-1];x<-as.matrix(dt[,-1]);rm(dt);calculate(path,genes,cells,x=x)
}else stop('Unsupported dataset')
stopifnot(length(scores)>0L)
ids<-unlist(lapply(scores,function(x)x$Source_Object_Cell_ID),use.names=FALSE)
stopifnot(!anyDuplicated(ids),setequal(ids,target))
fwrite(rbindlist(inventory),file.path(out,paste0(dataset,'_matrix_inventory.tsv')),sep='\t',na='',quote=FALSE)
p<-file.path(out,paste0(dataset,'_source_qc.tsv.gz'));tmp<-sub('.tsv.gz$','.partial.tsv.gz',p)
fwrite(rbindlist(scores),tmp,sep='\t',quote=FALSE,na='');stopifnot(file.rename(tmp,p))
message(Sys.time(),' SOURCE_QC_COMPLETE ',dataset)
