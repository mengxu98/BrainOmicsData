#!/usr/bin/env Rscript
# Full-cohort feature maps from observed counts and full-library denominators.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scattermore)})
if (!requireNamespace("scop", quietly=TRUE) ||
    as.character(utils::packageVersion("scop")) != "0.9.2")
  stop("Install the unified environment, including scop 0.9.2, before drawing S4")
if (!all(c("title.color", "title.face") %in% names(formals(scop::FeatureDimPlot))))
  stop("The installed scop package does not match the pinned plotting interface")
setDTthreads(4)
root<-normalizePath(commandArgs(TRUE)[1]);out<-normalizePath(commandArgs(TRUE)[2])
source('functions/utils.R')
source('functions/umap_style.R')
spec<-fread(file.path(out,'marker_order.tsv'));features<-spec$Gene
stopifnot(length(features)==33L,!anyDuplicated(features))
annotation_file<-Sys.getenv('BRAINOMICS_ANNOTATION_TABLE',unset=file.path('results','annotation','cluster_annotation.tsv'))
ann<-fread(annotation_file)
if ('CellType' %in% names(ann) && !'Working_CellType' %in% names(ann)) ann[,Working_CellType:=CellType]
ann[,Colour:=unname(brainomics_celltype_colors[Working_CellType])]
palette<-setNames(ann$Colour[!duplicated(ann$Working_CellType)],ann$Working_CellType[!duplicated(ann$Working_CellType)])
stopifnot(all(spec$Marker_Group %in% names(palette)))
find_one <- function(filename) {
 candidates<-list.files(root,recursive=TRUE,full.names=TRUE)
 matches<-candidates[basename(candidates)==filename]
 if(length(matches)!=1L) stop('Expected one ',filename,' under ',root,'; found ',length(matches))
 matches[[1L]]
}
metadata_file<-Sys.getenv('BRAINOMICS_FIG2_METADATA_FILE',unset='')
if(!nzchar(metadata_file)) metadata_file<-find_one('core_metadata_minimal.tsv.gz')
embedding_file<-Sys.getenv('BRAINOMICS_RPCA_UMAP_FILE',unset='')
if(!nzchar(embedding_file)) embedding_file<-find_one('embedding_umap.rpca.rds')
m<-fread(metadata_file,select=c('Cells','Cluster','Dataset'))
e<-readRDS(embedding_file)
stopifnot(nrow(m)==2602031L,!anyDuplicated(m$Cells),identical(rownames(e),m$Cells),all(is.finite(e)))
expression<-matrix(NA_real_,nrow(m),length(features),dimnames=list(NULL,features));seen<-logical(nrow(m))
counts_dir<-Sys.getenv('BRAINOMICS_MARKER_COUNTS_DIR',unset='')
if(!nzchar(counts_dir)) {
 candidates<-list.dirs(root,recursive=TRUE,full.names=TRUE)
 candidates<-candidates[basename(candidates)=='full_cell_evidence']
 if(length(candidates)!=1L) stop('Set BRAINOMICS_MARKER_COUNTS_DIR; found ',length(candidates),' candidate directories')
 counts_dir<-candidates[[1L]]
}
files<-list.files(counts_dir,pattern='_panel_counts.rds$',full.names=TRUE)
stopifnot(length(files)==22L)
for(f in files) {
 p<-readRDS(f);ix<-match(p$Cell,m$Cells);g<-match(features,p$Symbol)
 stopifnot(!anyNA(ix),!anyNA(g),!anyDuplicated(ix),!any(seen[ix]),length(unique(m$Dataset[ix]))==1L,
   identical(colnames(p$Counts),p$Cell),length(p$Library_Size)==length(ix),all(p$Library_Size>0))
 counts<-t(p$Counts[g,,drop=FALSE]);stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==floor(counts)))
 expression[ix,]<-log1p(counts/p$Library_Size*10000)
 seen[ix]<-TRUE
 rm(p,counts);gc(FALSE)
}
stopifnot(all(seen),all(is.finite(expression)),all(expression>=0))
marker_rows<-rbindlist(lapply(features,function(g) data.table(Cluster=m$Cluster,v=expression[,g])[,.(Gene=g,Cells=.N,MeanLog=mean(v),PctPositive=mean(v>0)),by=Cluster]))
evidence_file<-Sys.getenv('BRAINOMICS_MARKER_EVIDENCE_FILE',
  unset=file.path('results','annotation','full_cluster_marker_evidence.tsv'))
ref<-fread(evidence_file)
check<-merge(marker_rows,ref[,.(Cluster,Gene,Available_Cells,ReferenceMean=MeanLog,ReferencePct=PctPositive)],by=c('Cluster','Gene'))
stopifnot(nrow(check)==75*length(features),all(check$Cells==check$Available_Cells),max(abs(check$MeanLog-check$ReferenceMean))<1e-10,max(abs(check$PctPositive-check$ReferencePct))<1e-10)
summary<-rbindlist(lapply(features,function(g) data.table(Gene=g,Cells=nrow(m),Positive=sum(expression[,g]>0),Zero=sum(expression[,g]==0),Minimum=min(expression[,g]),Maximum=max(expression[,g]),Mean=mean(expression[,g]))))
fwrite(summary,file.path(out,'feature_expression_summary.tsv'),sep='\t')
title_colors<-setNames(unname(palette[spec$Marker_Group]),features)
spec[,Title_Color:=unname(title_colors[Gene])];fwrite(spec,file.path(out,'marker_order.tsv'),sep='\t')
normalized<-Matrix::Matrix(t(expression),sparse=TRUE)
rownames(normalized)<-features;colnames(normalized)<-m$Cells
obj<-SeuratObject::CreateSeuratObject(counts=SeuratObject::CreateAssayObject(data=normalized),assay='RNA')
colnames(e)<-c('UMAP_1','UMAP_2')
obj[['umap']]<-SeuratObject::CreateDimReducObject(embeddings=e,key='UMAP_',assay='RNA')
stopifnot(identical(colnames(obj),m$Cells))
rm(normalized);gc(FALSE)
# Keep the legend compact while leaving readable labels and explicit gaps.
legend_scale <- 0.65
raster_pixels <- 600
raster_point_radius <- 0.6
figure_width_mm <- 240
figure_height_mm <- 270
plots<-vector('list',length(features))
cairo_pdf(file.path(out,'measurement.pdf'),width=8,height=6,family='Arial')
for(j in seq_along(features)) {
 gene<-features[j];v<-expression[,gene]
 q<-scop::FeatureDimPlot(obj,features=gene,reduction='umap',layer='data',
   show_stat=FALSE,palette='Spectral',bg_cutoff=0,bg_color='grey80',
   lower_quantile=0,upper_quantile=1,lower_cutoff=0,keep_scale='feature',
   raster=TRUE,raster.dpi=rep(raster_pixels,2),pt.size=raster_point_radius*512/raster_pixels,
   title.color=title_colors,title.face=NULL,subtitle=NULL,
   legend.position='bottom',legend.direction='horizontal',legend.title='log1p(CP10K)',
   theme_use=theme_blank_axis,theme_args=modifyList(umap_axis_args,list(add_coord=FALSE,xlab='',ylab='',lab_size=8.5,axis_lwd=.5)),
   combine=FALSE,seed=20260918,verbose=FALSE)[[1]]+
   theme(text=element_text(family='Arial'),strip.text.x=element_text(size=10.5),
     strip.background=element_blank(),legend.title=element_text(size=8.5*legend_scale,margin=margin(b=1.5,unit='pt')),
     legend.text=element_text(size=8*legend_scale,margin=margin(t=1.5,unit='pt')),
     legend.margin=margin(0,0,0,0,'mm'),
     legend.box.spacing=grid::unit(0.5,'mm'),plot.margin=margin(0.25,0.25,0.5,0.25,'mm'))+
   guides(colour=guide_colourbar(title.position='top',barwidth=grid::unit(26*legend_scale,'mm'),barheight=grid::unit(2.5*legend_scale,'mm'),frame.colour='black',ticks.colour='black',
     theme=theme(legend.text=element_text(size=8*legend_scale,margin=margin(t=1.5,unit='pt')),
       legend.title=element_text(size=8.5*legend_scale,margin=margin(b=1.5,unit='pt')),
       legend.ticks.length=grid::unit(0.25,'mm'),legend.frame=element_rect(linewidth=0.25),legend.ticks=element_line(linewidth=0.2))))
 ix<-match(rownames(q$data),m$Cells)
 stopifnot(nrow(q$data)==nrow(m),!anyNA(ix),!anyDuplicated(ix),
   all(is.na(q$data$value)==(v[ix]==0)),all(q$data$value[!is.na(q$data$value)]==v[ix][v[ix]>0]),
   identical(q$theme$strip.text.x$colour,unname(title_colors[gene])),identical(q$theme$strip.text.x$face,'italic'),
   is.null(q$labels$subtitle))
 plots[[j]]<-ggplotGrob(q);
 rm(q);gc(FALSE)
}
shared_heights<-do.call(grid::unit.pmax,lapply(plots,function(g)g$heights));plots<-lapply(plots,function(g){g$heights<-shared_heights;g})
p<-wrap_plots(c(lapply(plots,wrap_elements),rep(list(plot_spacer()),3)),ncol=6) &
  theme(plot.margin=margin(0,0,0,0,'mm'))
g<-patchworkGrob(p)
cairo_pdf(file.path(out,'figS4.pdf'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,family='Arial');grid::grid.draw(g);dev.off()
dev.off();unlink(file.path(out,'measurement.pdf'))
message('Supplementary Figure S4 written')
