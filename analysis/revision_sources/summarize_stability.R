suppressPackageStartupMessages({library(jsonlite)})
a<-commandArgs(TRUE);stopifnot(length(a)==2L);base<-normalizePath(a[1]);out<-normalizePath(a[2])
old<-readRDS(file.path(base,'integration_preview_20260917/cluster_assignments.rds'))
formal<-read.delim(file.path(out,'final_cluster_annotation.tsv'),stringsAsFactors=FALSE)
stopifnot(nrow(old)==2602031L,nrow(formal)==75L,length(unique(formal$CellType))==12L)
dataset<-sub('::.*$','',old$Cell);stopifnot(length(unique(dataset))==22L)
class<-formal$CellType[match(old$Cluster,formal$Cluster)];stopifnot(!anyNA(class))
ari<-function(x,y){t<-table(x,y);c2<-function(z) z*(z-1)/2;n<-sum(t);a<-sum(c2(rowSums(t)));b<-sum(c2(colSums(t)));e<-a*b/c2(n);den<-(a+b)/2-e;if(den==0)return(if(identical(x,y))1 else NA_real_);(sum(c2(t))-e)/den}
res<-list();studies<-list();classes<-list();clusters<-list()
for(i in 1:5){
 run<-file.path(out,paste0('run_',i));stopifnot(file.exists(file.path(run,'COMPLETE.json')))
 z<-readRDS(file.path(run,'assignments.rds'));stopifnot(identical(z$Cell,old$Cell))
 s<-read.delim(file.path(run,'summary.tsv'));res[[i]]<-s
 counts<-table(class,z$Cluster);majority<-rownames(counts)[max.col(t(counts),ties.method='first')];names(majority)<-colnames(counts);mapped<-unname(majority[z$Cluster])
 classes[[i]]<-do.call(rbind,lapply(sort(unique(class)),function(k){ix<-class==k;data.frame(Task=i,CellType=k,Cells=sum(ix),Majority_Mapped_Agreement=mean(mapped[ix]==class[ix]))}))
 studies[[i]]<-do.call(rbind,lapply(sort(unique(dataset)),function(k){ix<-dataset==k;data.frame(Task=i,Dataset=k,Cells=sum(ix),Partition_ARI=ari(old$Cluster[ix],z$Cluster[ix]),Majority_Mapped_Class_Agreement=mean(mapped[ix]==class[ix]))}))
 c<-read.delim(file.path(run,'cluster_stability.tsv'));c$Task<-i;clusters[[i]]<-c
 rm(z,mapped);gc()
}
for(x in list(list(name='summary',value=res),list(name='by_dataset',value=studies),list(name='by_celltype',value=classes),list(name='by_cluster',value=clusters)))write.table(do.call(rbind,x$value),file.path(out,paste0('stability_',x$name,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(state='PASS',runs=5,cells=2602031,datasets=22,baseline_replay=TRUE,scope='All cells on one fixed production graph. Majority label matching measures stability conditional on the adopted annotation and is not independent biological accuracy. Integration randomness and graph construction were not perturbed.'),file.path(out,'SUMMARY_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
