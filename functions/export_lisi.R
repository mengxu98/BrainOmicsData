#!/usr/bin/env Rscript
# Serialize completed per-cell LISI scores; no neighbor search or recomputation.
#
# Argument 4 (optional) names the dataset-LISI result directory.  Both layouts are
# accepted: <dir>/<Method>/<reduction>_dataset_lisi.rds and the flat
# <dir>/<method>_dataset_lisi.rds.  A score object may be a list with Cell/Score or
# a named numeric vector (names are the cell identifiers).  Only the cell order and
# score values are validated; existing results are never recomputed.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(4L)
a<-commandArgs(TRUE);stopifnot(length(a)%in%c(3L,4L))
root<-normalizePath(a[1]);package<-a[2];internal<-a[3]
# An existing per-cell LISI result is authoritative: keep it and stop.  Result
# files are only rebuilt when none is present and a checkpoint directory is given.
existing<-file.path(package,'validation/lisi.tsv.gz')
if(file.exists(existing)){
  check<-fread(cmd=paste('gzip -cd',shQuote(existing)),sep='\t',nrows=2L)
  if(ncol(check)==9L){
    rows<-as.integer(system2('bash',c('-c',shQuote(paste('gzip -cd',shQuote(existing),'| wc -l'))),stdout=TRUE))
    if(!is.na(rows)&&rows-1L==2602031L){
      write_json(list(state='KEPT_EXISTING_RESULT',cells=2602031L,columns=names(check),recomputed=FALSE),
                 file.path(internal,'lisi_export.json'),pretty=TRUE,auto_unbox=TRUE)
      message(Sys.time(),' existing validation/lisi.tsv.gz kept; nothing recomputed')
      quit(save='no',status=0L)
    }
  }
  stop('Existing validation/lisi.tsv.gz does not match the expected 2,602,031 x 9 layout')
}
frozen<-normalizePath(Sys.getenv('BRAINOMICS_FROZEN_DIR',unset=file.path(dirname(root),'frozen_run')),mustWork=TRUE)
cl<-readRDS(file.path(frozen,'cluster_assignments.rds'))
stopifnot(nrow(cl)==2602031L,!anyDuplicated(cl$Cell))
score_source<-function(dir,method,reduction){
  candidates<-c(file.path(dir,method,paste0(reduction,'_dataset_lisi.rds')),
                file.path(dir,paste0(tolower(method),'_dataset_lisi.rds')),
                file.path(dir,paste0(reduction,'_dataset_lisi.rds')))
  hit<-candidates[file.exists(candidates)]
  if(length(hit)==0L) stop('No dataset-LISI result for ',method,' under ',dir)
  hit[[1L]]
}
score_values<-function(path){
  score<-readRDS(path)
  if(is.list(score)&&all(c('Cell','Score')%in%names(score))) return(list(Cell=score$Cell,Score=score$Score))
  if(is.numeric(score)&&!is.null(names(score))) return(list(Cell=names(score),Score=unname(score)))
  stop('Unsupported score object: ',path)
}
v<-data.table(Cells=sprintf('Cell%07d',seq_len(nrow(cl))))
methods<-c(Raw='pca',scVI='integrated.scvi',Harmony='integrated.harmony',RPCA='integrated.rpca')
summaries<-list()
for(method in names(methods)){
  message(Sys.time(),' export saved scores ',method)
  dsfile<-if(length(a)==4L) score_source(normalizePath(a[4],mustWork=TRUE),method,methods[[method]]) else stop('Supply the dataset-LISI result directory as argument 4')
  ds<-score_values(dsfile)
  at<-match(cl$Cell,as.character(ds$Cell))
  if(anyNA(at)||anyDuplicated(ds$Cell)) stop('Score cells differ from the released cohort: ',dsfile)
  ds$Score<-ds$Score[at]
  csfile<-file.path(root,'07_downstream/revision_20260918/25_formal_manuscript_review_20260918/tables/full_lisi/latent50',method,'cell_scores.rds')
  cs<-readRDS(csfile)
  stopifnot(identical(as.character(cs$Cells),as.character(cl$Cell)),
            length(ds$Score)==nrow(cl),nrow(cs$cLISI)==nrow(cl),all(is.finite(ds$Score)),all(is.finite(cs$cLISI[,'Source_Full'])))
  name<-if(method=='Raw')'Raw_PCA' else method
  v[,(paste0(name,'_iLISI')):=ds$Score]
  v[,(paste0(name,'_cLISI')):=cs$cLISI[,'Source_Full']]
  summaries[[method]]<-list(method=method,cells=nrow(v),dataset_scores_source=dsfile,source_annotation_scores_source=csfile,cell_order_identical=TRUE)
  rm(ds,cs);gc()
}
setcolorder(v,c('Cells',paste0(c('Raw_PCA','scVI','Harmony','RPCA'),'_iLISI'),paste0(c('Raw_PCA','scVI','Harmony','RPCA'),'_cLISI')))
dir.create(file.path(package,'validation'),showWarnings=FALSE)
tmp<-file.path(package,'validation/lisi.partial');fwrite(v,tmp,sep='\t',quote=FALSE,na='',nThread=4)
stopifnot(system2('pigz',c('-n','-f','-p','4',shQuote(tmp)))==0L,file.rename(paste0(tmp,'.gz'),file.path(package,'validation/lisi.tsv.gz')))
write_json(list(state='EXPORTED',cells=nrow(v),columns=names(v),completed_score_inputs=summaries,recomputed=FALSE),file.path(internal,'lisi_export.json'),pretty=TRUE,auto_unbox=TRUE)
message(Sys.time(),' LISI_EXPORT_COMPLETE')
