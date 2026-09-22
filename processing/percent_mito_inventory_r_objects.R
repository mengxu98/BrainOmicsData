#!/usr/bin/env Rscript
# Count mitochondrial reads in the original source objects (RDS/RData).
# Source-cell identifiers are kept as they appear; map_to_cohort.R joins them.
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(data.table)})
setDTthreads(2L)
a<-commandArgs(TRUE);stopifnot(length(a)==3L);raw<-a[1];out<-a[2];dataset<-a[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files<-list(GSE204683='GSE204683/GSE204683_count_matrix.RDS',HYPOMAP='HYPOMAP/human_HYPOMAP_snRNASeq.rds',Ma_et_al_2022='GSE207334/Ma_Sestan_mat.rds',SomaMut='SomaMut/pfc.clean.rds',GSE296073='GSE296073/h_pre_peri_DY.rds',Li_et_al_2018=c('Li_et_al_2018/rawData/Prenatal_scRNA_seq/Sestan.fetalHuman.Psychencode.Rdata','Li_et_al_2018/rawData/Adult_snRNA_seq/Sestan.adultHumanNuclei.Psychencode.Rdata'))
stopifnot(dataset%in%names(files));inventory<-list();scores<-list()
inspect<-function(x,path,layer){
 stopifnot(length(dim(x))==2L,!is.null(rownames(x)),!is.null(colnames(x)),!anyDuplicated(colnames(x)))
 genes<-sub('^.*[|]','',rownames(x));mt<-grepl('^MT-',genes);nmt<-sum(mt)
 state<-if(nmt>0L)'COMPUTED_SOURCE_COUNTS' else if(any(grepl('^ENSG[0-9]',genes))) 'UNRESOLVED_GENE_IDS' else 'NO_MT_FEATURES'
 message(Sys.time(),' ',dataset,' ',layer,' ',nrow(x),' features; MT=',nmt,'; ',ncol(x),' cells')
 if(nmt>0L){
  if(inherits(x,'sparseMatrix')){
   xx<-if('x'%in%slotNames(x)) x@x else stop('Non-numeric sparse matrix')
   if(length(xx))for(i in seq.int(1L,length(xx),by=10000000L)){v<-xx[i:min(length(xx),i+9999999L)];stopifnot(all(is.finite(v)),all(v>=0),all(v==floor(v)))}
  }else for(i in seq_len(ncol(x))){v<-x[,i];stopifnot(all(is.finite(v)),all(v>=0),all(v==floor(v)))}
  total<-Matrix::colSums(x);numerator<-Matrix::colSums(x[mt,,drop=FALSE]);pct<-ifelse(total>0,100*numerator/total,NA_real_)
  stopifnot(all(is.na(pct)|(is.finite(pct)&pct>=0&pct<=100)))
 }else{total<-numerator<-pct<-rep(NA_real_,ncol(x))}
 key<-paste0(basename(path),'::',layer)
 scores[[key]]<<-data.table(Dataset=dataset,Source_Object_Cell_ID=colnames(x),percent_mito=pct,MT_Counts=numerator,Total_Counts=total,Source_File=path,Source_Layer=layer,Status=if(nmt>0L)ifelse(is.na(pct),'ZERO_TOTAL',state)else state)
 inventory[[key]]<<-data.table(Dataset=dataset,Source_File=path,Kind='R_object',Source_Bytes=file.info(path)$size,Count_Semantics='Original source count slot; exon/intron definition retained from source, not inferred',Features=nrow(x),MT_Genes=nmt,MT_Symbols=paste(sort(unique(genes[mt])),collapse=';'),Status=state,Source_Layer=layer,Source_Cells=ncol(x))
}
for(rel in files[[dataset]]){
 path<-file.path(raw,rel);message(Sys.time(),' Loading ',path)
 if(grepl('Rdata$',path)){
  env<-new.env();load(path,envir=env);field<-if(grepl('fetalHuman',path))'count2' else 'umi.raw';inspect(env[[field]],path,field);rm(env)
 }else{
  obj<-readRDS(path)
  if(inherits(obj,'Seurat')){
   assay<-if('RNA'%in%names(obj@assays))'RNA' else stop('No explicit RNA assay; inspect manually')
   if(inherits(obj@assays[[assay]],'Assay5'))for(layer in Layers(obj@assays[[assay]],search='^counts'))inspect(LayerData(obj,assay=assay,layer=layer),path,paste(assay,layer,sep='/')) else inspect(obj@assays[[assay]]@counts,path,'RNA/counts')
  }else inspect(obj,path,'matrix')
  rm(obj)
 }
 gc()
}
stopifnot(length(inventory)>0L)
fwrite(rbindlist(inventory,fill=TRUE),file.path(out,paste0(dataset,'_inventory.tsv')),sep='\t',quote=FALSE,na='')
p<-file.path(out,paste0(dataset,'_source_qc.tsv.gz'));tmp<-sub('.tsv.gz$','.partial.tsv.gz',p)
fwrite(rbindlist(scores),tmp,sep='\t',quote=FALSE,na='');stopifnot(file.rename(tmp,p))
message(Sys.time(),' SOURCE_QC_COMPLETE ',dataset)
