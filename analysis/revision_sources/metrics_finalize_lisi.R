suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest)})
p<-commandArgs(TRUE)[1];out<-file.path(p,'tables/full_lisi');dsall<-list();eff<-list();cells<-list();setDTthreads(4)
for(space in c('latent50','umap2'))for(scheme in c('Dataset','Source_Full','Formal_Type')){
 dd<-lapply(c('Raw','RPCA','Harmony','scVI'),function(method){path<-file.path(out,space,method);stopifnot(file.exists(file.path(path,'COMPLETE.json')));x<-fread(file.path(path,paste0(scheme,'_donor_summary.tsv')));x[,Method:=method];x});names(dd)<-c('Raw','RPCA','Harmony','scVI')
 for(method in names(dd)){
 x<-merge(dd[[method]],dd$Raw[,.(Dataset,Canonical_Donor_ID,Raw_Score=cLISI)],by=c('Dataset','Canonical_Donor_ID'));x[,Change_From_Raw:=cLISI-Raw_Score];ds<-x[,.(Cells=sum(Cells),Donors=.N,Mean=mean(cLISI),Change_From_Raw=mean(Change_From_Raw),Neighbor30_Purity=mean(Neighbor30_Purity)),by=Dataset];ds[,`:=`(Method=method,Space=space,Label_Scheme=scheme)];dsall[[paste(space,scheme,method)]]<-ds
 for(sens in c('all_studies','exclude_Ma','exclude_three_linked')){
 z<-switch(sens,all_studies=ds,exclude_Ma=ds[Dataset!='Ma_et_al_2022'],exclude_three_linked=ds[!Dataset%in%c('Li_et_al_2018','Ma_et_al_2022','GSE186538')]);h<-qt(.975,nrow(z)-1)*sd(z$Change_From_Raw)/sqrt(nrow(z));eff[[paste(space,scheme,method,sens)]]<-data.table(Space=space,Label_Scheme=scheme,Method=method,Sensitivity=sens,Studies=nrow(z),Study_Equal_Mean=mean(z$Mean),Change_From_Raw=mean(z$Change_From_Raw),CI95_Lower=mean(z$Change_From_Raw)-h,CI95_Upper=mean(z$Change_From_Raw)+h)
 }
 }
}
for(space in c('latent50','umap2'))for(method in c('Raw','RPCA','Harmony','scVI')){x<-fread(file.path(out,space,method,'summary.tsv'));x[,Space:=space];cells[[paste(space,method)]]<-x}
fwrite(rbindlist(dsall),file.path(p,'tables/full_lisi_dataset_summary.tsv'),sep='\t');fwrite(rbindlist(eff),file.path(p,'tables/full_lisi_study_effects.tsv'),sep='\t');fwrite(rbindlist(cells),file.path(p,'tables/full_lisi_cell_summary.tsv'),sep='\t');write_json(list(state='COMPLETE',methods=4,spaces=2,cells_each=2602031,schemes=c('Dataset','Source_Full','Formal_Type'),sampling=FALSE,resampling=FALSE,aggregation='Cell -> canonical donor within study -> equal study; paired changes versus Raw; pointwise Student t CI across study means',HNSW=list(M=24,construction_ef=200,search_ef=400,seed=2026),perplexity=30,neighbors=89,limits='Approximate HNSW, not RANN equivalence; source categories include unresolved and cycling labels; Formal_Type depends on RPCA annotation'),file.path(out,'COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
