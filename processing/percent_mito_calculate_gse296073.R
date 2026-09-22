#!/usr/bin/env Rscript
# Mitochondrial fraction for GSE296073 from its 28 original 10x libraries.
# The author RDS for this source carries no MT genes.
#
# Usage: percent_mito_calculate_gse296073.R RAW_ROOT OUTPUT_DIR GSE296073
suppressPackageStartupMessages({library(Matrix);library(data.table)})
setDTthreads(2L)
a<-commandArgs(TRUE);stopifnot(length(a)==3L);raw<-a[1];out<-a[2];dataset<-a[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE);scores<-list();inventory<-list()
calculate<-function(path,genes,cells,keep=seq_along(genes),x=NULL){
 message(Sys.time(),' Reading ',path)
 if(is.null(x)){con<-gzfile(path,'rt');x<-readMM(con);close(con)}
 stopifnot(nrow(x)==length(genes),ncol(x)==length(cells),!anyDuplicated(cells))
 x<-x[keep,,drop=FALSE];genes<-genes[keep];mt<-grepl('^MT-',sub('^.*[|]','',genes));stopifnot(any(mt))
 if(inherits(x,'sparseMatrix')){
  values<-x@x
  if(length(values))for(i in seq.int(1L,length(values),by=10000000L)){v<-values[i:min(length(values),i+9999999L)];stopifnot(all(is.finite(v)),all(v>=0),all(v==floor(v)))}
 }else for(i in seq_len(ncol(x))){v<-x[,i];stopifnot(all(is.finite(v)),all(v>=0),all(v==floor(v)))}
 total<-Matrix::colSums(x);numerator<-Matrix::colSums(x[mt,,drop=FALSE]);pct<-ifelse(total>0,100*numerator/total,NA_real_)
 stopifnot(all(is.na(pct)|(is.finite(pct)&pct>=0&pct<=100)))
 scores[[path]]<<-data.table(Dataset=dataset,Source_Object_Cell_ID=cells,percent_mito=pct,MT_Counts=numerator,Total_Counts=total,Source_File=path,Source_Layer='source Gene Expression counts',Status=ifelse(total>0,'COMPUTED_SOURCE_COUNTS','ZERO_TOTAL'))
 inventory[[path]]<<-data.table(Dataset=dataset,Source_File=path,Kind='source_count_matrix',Source_Bytes=file.info(path)$size,Count_Semantics='Unmodified source Gene Expression feature universe; exon/intron scope as supplied, not independently verified',Features=length(genes),MT_Genes=sum(mt),MT_Symbols=paste(sort(unique(genes[mt])),collapse=';'),Status='COMPUTED_SOURCE_COUNTS',Source_Cells=length(cells))
 message(Sys.time(),' Computed ',length(cells),' source cells; MT=',sum(mt));rm(x);gc()
}
stopifnot(dataset=='GSE296073')
meta<-fread(file.path(dirname(raw),'processed/GSE296073/metadata_raw.tsv.gz'),select=c('Source_Object_Cell_ID','libraryID','age'))
libs<-unique(meta[,.(libraryID,age)])
libs[,Source_Suffix:=mapply(function(x,a)sub(paste0(a,'_'),'_',sub('^h','',x)),libraryID,age)]
stopifnot(nrow(libs)==28L,!anyDuplicated(libs$Source_Suffix))
fs<-list.files(file.path(raw,dataset),'features.tsv.gz$',recursive=TRUE,full.names=TRUE)
suffix<-sub('^GSM[0-9]+_','',basename(dirname(fs)))
stopifnot(length(fs)==28L,setequal(suffix,libs$Source_Suffix))
for(k in seq_along(fs)){
 f<-fs[k];lib<-libs$libraryID[match(suffix[k],libs$Source_Suffix)]
 feat<-fread(f,header=FALSE);genes<-feat[[2]];keep<-if(ncol(feat)>=3L)which(feat[[3]]=='Gene Expression')else seq_along(genes)
 cells<-paste0(lib,'_',fread(sub('features.tsv.gz$','barcodes.tsv.gz',f),header=FALSE)[[1]])
 calculate(sub('features.tsv.gz$','matrix.mtx.gz',f),genes,cells,keep)
}
all_ids<-unlist(lapply(scores,function(x)x$Source_Object_Cell_ID),use.names=FALSE)
stopifnot(!anyDuplicated(all_ids))
missing<-meta[!Source_Object_Cell_ID%in%all_ids]
fwrite(missing,file.path(out,'GSE296073_cells_absent_from_10x.tsv.gz'),sep='\t',quote=FALSE)
if(nrow(missing))scores[['source_cells_absent']]<-data.table(Dataset=dataset,Source_Object_Cell_ID=missing$Source_Object_Cell_ID,percent_mito=NA_real_,MT_Counts=NA_real_,Total_Counts=NA_real_,Source_File=file.path(raw,dataset),Source_Layer='28 original 10x libraries; exact library+barcode match only',Status='SOURCE_CELL_ABSENT')
message('Exact source-library join: ',length(all_ids),' available barcodes; ',nrow(missing),' source-object cells absent. No cross-library matching or zero imputation.')
fwrite(libs,file.path(out,'GSE296073_library_crosswalk.tsv'),sep='\t',quote=FALSE)
stopifnot(length(scores)>0L)
fwrite(rbindlist(inventory),file.path(out,paste0(dataset,'_matrix_inventory.tsv')),sep='\t',na='',quote=FALSE)
p<-file.path(out,paste0(dataset,'_source_qc.tsv.gz'));tmp<-sub('.tsv.gz$','.partial.tsv.gz',p)
fwrite(rbindlist(scores),tmp,sep='\t',quote=FALSE,na='');stopifnot(file.rename(tmp,p))
message(Sys.time(),' SOURCE_QC_COMPLETE ',dataset)
