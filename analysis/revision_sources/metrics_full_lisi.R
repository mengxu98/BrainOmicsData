#!/usr/bin/env Rscript
args<-commandArgs(TRUE);root<-normalizePath(args[1]);scope<-normalizePath(args[2]);method<-args[3];space<-args[4]
run<-file.path(root,'07_downstream/revision_20260918');out<-file.path(scope,'tables/full_lisi',space,method);dir.create(out,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(run,'environment/R-library'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(RcppHNSW);library(lisi);library(jsonlite);library(Rcpp);library(digest)})
setDTthreads(4)
pin<-readLines(file.path(find.package('lisi'),'BRAINOMICS_PINNED_SOURCE'));stopifnot('numeric_patch=double_precision_probability_normalization_v1'%in%pin,'numeric_patch_sha256=6c8d1e25cf62aaa385ef9fff90eeedfb31519768a9853c5d5a116a85837171dd'%in%pin)
writeLines(pin,file.path(out,'lisi_source_pin.txt'))
message(Sys.time(),' loading full metadata ',method)
m<-as.data.table(readRDS(file.path(run,'01_metadata/metadata_working.rds')))
lab<-as.data.table(readRDS(file.path(root,'01_integration_evaluation/independent_labels/independent_labels_by_cell.rds')))
j<-match(m$Cells,lab$Cells);stopifnot(!anyNA(j),!anyDuplicated(m$Cells),nrow(m)==2602031L);lab<-lab[j]
m[,Source_Full:=as.character(lab$Source_Coarse)]
# Retain every cell; restore the author-provided TAC cycling state, without inventing a lineage.
tac<-which(m$Dataset=='GSE217511' & grepl('^TAC',m$CellType_raw));stopifnot(length(tac)==1209L,all(m$Source_Full[tac]=='Unassigned'))
m[tac,Source_Full:='TAC_cycling_state'];m[is.na(Source_Full)|Source_Full=='Unassigned',Source_Full:='Unresolved_source_label']
stopifnot(sum(m$Source_Full=='Unresolved_source_label')==31L)
a<-fread(file.path(run,'18_final_annotation_20260918/cluster_annotation.tsv'));m[,Formal_Type:=a$Working_CellType[match(Cluster,a$Cluster)]]
schemes<-c('Source_Full','Formal_Type','Dataset');labels<-lapply(schemes,function(s)as.integer(factor(m[[s]]))-1L);names(labels)<-schemes;levels_n<-sapply(labels,function(x)length(unique(x)));rm(lab,a);gc(FALSE)
cppFunction('List drop_self(IntegerMatrix ix, NumericMatrix ds, int start) { int n=ix.nrow(); IntegerMatrix oi(n,89); NumericMatrix od(n,89); for(int i=0;i<n;i++){int q=0;for(int j=0;j<ix.ncol() && q<89;j++){if(ix(i,j)==start+i)continue;oi(i,q)=ix(i,j);od(i,q)=ds(i,j);q++;}if(q!=89)stop("insufficient nonself neighbors");}return List::create(Named("idx")=oi,Named("dist")=od);}')
methods<-c(Raw='pca',RPCA='integrated.rpca',Harmony='integrated.harmony',scVI='integrated.scvi');stopifnot(method%in%names(methods))
if(space=='umap2')methods<-c(Raw='umap.unintegrated',RPCA='umap.rpca',Harmony='umap.harmony',scVI='umap.scvi')
efile<-file.path(root,'00_input_audit/compact',paste0('embedding_',methods[[method]],'.rds'));X<-readRDS(efile);stopifnot(identical(rownames(X),m$Cells),all(is.finite(X)))
message(Sys.time(),' building full-cohort graph ',method,'; cells=',nrow(X))
index<-hnsw_build(X,distance='euclidean',M=24,ef=200,n_threads=4,random_seed=2026)
value<-matrix(NA_real_,nrow(m),length(schemes),dimnames=list(NULL,schemes));purity<-value
for(st in seq.int(1L,nrow(m),by=50000L)){
 en<-min(st+49999L,nrow(m));ii<-st:en;h<-hnsw_search(X[ii,,drop=FALSE],index,k=90L,ef=400,n_threads=4);h<-drop_self(h$idx,h$dist,st)
 ids0<-t(h$idx)-1L;ds<-t(h$dist)
 for(s in schemes){
  v<-1/lisi:::compute_simpson_index(ds,ids0,labels[[s]],levels_n[[s]],30)
  stopifnot(all(is.finite(v)),all(v>=1-1e-12),all(v<=levels_n[[s]]+1e-12));value[ii,s]<-pmax(v,1)
  purity[ii,s]<-rowMeans(matrix(labels[[s]][h$idx[,1:30,drop=FALSE]],nrow=length(ii))==labels[[s]][ii])
 }
 message(Sys.time(),' ',method,' scored ',en,'/',nrow(m),' all 3 label schemes')
 write_json(list(method=method,state='RUNNING',scored_cells=en,total_cells=nrow(m),updated=as.character(Sys.time())),file.path(out,'progress.json'),auto_unbox=TRUE,pretty=TRUE)
}
stopifnot(!anyNA(value),nrow(value)==2602031L)
saveRDS(list(Cells=m$Cells,cLISI=value,Neighbor30_Purity=purity,method=method,sampling=FALSE),file.path(out,'cell_scores.rds'))
res<-list()
for(s in schemes){
 m[,`:=`(cLISI=value[,s],Neighbor30_Purity=purity[,s])]
 d<-m[,.(Cells=.N,cLISI=mean(cLISI),Neighbor30_Purity=mean(Neighbor30_Purity)),by=.(Dataset,Canonical_Donor_ID)]
 fwrite(d,file.path(out,paste0(s,'_donor_summary.tsv')),sep='\t')
 fwrite(d[,.(Donors=.N,Cells=sum(Cells),cLISI=mean(cLISI),Neighbor30_Purity=mean(Neighbor30_Purity)),by=Dataset],file.path(out,paste0(s,'_dataset_summary.tsv')),sep='\t')
 fwrite(m[,.(Cells=.N,cLISI=mean(cLISI),Neighbor30_Purity=mean(Neighbor30_Purity)),by=c('Cluster',s)],file.path(out,paste0(s,'_cluster_summary.tsv')),sep='\t')
 res[[s]]<-data.table(Method=method,Label_Scheme=s,Cells=nrow(m),Label_Categories=levels_n[[s]],Mean=mean(value[,s]),Median=median(value[,s]),Q25=unname(quantile(value[,s],.25)),Q75=unname(quantile(value[,s],.75)))
}
fwrite(rbindlist(res),file.path(out,'summary.tsv'),sep='\t')
write_json(list(state='COMPLETE',method=method,space=space,formal_annotation_sha256=digest(file=file.path(run,'18_final_annotation_20260918/cluster_annotation.tsv'),algo='sha256'),input_cells=nrow(m),query_cells=nrow(m),neighbor_pool_cells=nrow(X),sampling=FALSE,excluded_cells=0,perplexity=30,neighbors_excluding_self=89,HNSW=list(M=24,construction_ef=200,search_ef=400,threads=4,seed=2026),label_schemes=schemes,source_unresolved_cells=31,source_TAC_cycling_state_cells=1209,embedding_sha256=digest(file=efile,algo='sha256'),script_sha256=digest(file=file.path(scope,'scripts/metrics_full_lisi.R'),algo='sha256'),limitations='HNSW approximate neighbors, no exact-recall validation claimed. Source labels include differing granularity, a cycling state and 31 unresolved labels. Working/candidate cLISI depends on RPCA-derived annotation and is not independent validation of RPCA.'),file.path(out,'COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
message(Sys.time(),' COMPLETE ',method)
