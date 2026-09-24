#!/usr/bin/env Rscript
# Full-cohort feature maps from observed counts and full-library denominators.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scattermore);library(jsonlite);library(digest)})
scop_source <- Sys.getenv("SCOP_SOURCE_PATH")
if (nzchar(scop_source)) {
  pkgload::load_all(scop_source, quiet=TRUE, compile=FALSE)
} else if (!requireNamespace("scop", quietly=TRUE)) {
  stop("Install scop or set SCOP_SOURCE_PATH to its checkout")
}
if (!all(c("title.color", "title.face") %in% names(formals(scop::FeatureDimPlot))))
  stop("This scop build lacks title.color/title.face; set SCOP_SOURCE_PATH to the updated checkout")
setDTthreads(4)
root<-normalizePath(commandArgs(TRUE)[1]);out<-normalizePath(commandArgs(TRUE)[2])
source('functions/utils.R')
source('functions/umap_style.R')
spec<-fread(file.path(out,'marker_order.tsv'));features<-spec$Gene
stopifnot(length(features)==33L,!anyDuplicated(features))
ann<-fread(file.path(root,'07_downstream/revision_20260918/11_c65_microglia_compact_20260918/cluster_annotation.tsv'))
palette<-setNames(ann$Colour[!duplicated(ann$Working_CellType)],ann$Working_CellType[!duplicated(ann$Working_CellType)])
stopifnot(all(spec$Marker_Group %in% names(palette)))
m<-fread(file.path(root,'00_input_audit/compact/core_metadata_minimal.tsv.gz'),select=c('Cells','Cluster','Dataset'))
e<-readRDS(file.path(root,'00_input_audit/compact/embedding_umap.rpca.rds'))
stopifnot(nrow(m)==2602031L,!anyDuplicated(m$Cells),identical(rownames(e),m$Cells),all(is.finite(e)))
expression<-matrix(NA_real_,nrow(m),length(features),dimnames=list(NULL,features));seen<-logical(nrow(m));inputs<-list()
files<-list.files(file.path(root,'03_annotation/formal_20260917/full_cell_evidence'),pattern='_panel_counts.rds$',full.names=TRUE)
stopifnot(length(files)==22L)
for(f in files) {
 p<-readRDS(f);ix<-match(p$Cell,m$Cells);g<-match(features,p$Symbol)
 stopifnot(!anyNA(ix),!anyNA(g),!anyDuplicated(ix),!any(seen[ix]),length(unique(m$Dataset[ix]))==1L,
   identical(colnames(p$Counts),p$Cell),length(p$Library_Size)==length(ix),all(p$Library_Size>0))
 counts<-t(p$Counts[g,,drop=FALSE]);stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==floor(counts)))
 expression[ix,]<-log1p(counts/p$Library_Size*10000)
 seen[ix]<-TRUE
 inputs[[length(inputs)+1L]]<-data.table(File=basename(f),Dataset=unique(m$Dataset[ix]),Cells=length(ix),SHA256=digest(file=f,algo='sha256'))
 message(Sys.time(),' loaded ',basename(f),' ',length(ix),' cells')
 rm(p,counts);gc(FALSE)
}
stopifnot(all(seen),all(is.finite(expression)),all(expression>=0))
fwrite(rbindlist(inputs),file.path(out,'input_manifest.tsv'),sep='\t')
# Verify the reconstructed expression against the existing complete-cluster evidence.
marker_rows<-rbindlist(lapply(features,function(g) data.table(Cluster=m$Cluster,v=expression[,g])[,.(Gene=g,Cells=.N,MeanLog=mean(v),PctPositive=mean(v>0)),by=Cluster]))
ref<-fread(file.path(root,'07_downstream/revision_20260918/11_c65_microglia_compact_20260918/outputs/annotation/full_cluster_marker_evidence.tsv'))
check<-merge(marker_rows,ref[,.(Cluster,Gene,Available_Cells,ReferenceMean=MeanLog,ReferencePct=PctPositive)],by=c('Cluster','Gene'))
stopifnot(nrow(check)==75*length(features),all(check$Cells==check$Available_Cells),max(abs(check$MeanLog-check$ReferenceMean))<1e-10,max(abs(check$PctPositive-check$ReferencePct))<1e-10)
fwrite(check,file.path(out,'full_cluster_expression_verification.tsv'),sep='\t')
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
plots<-vector('list',length(features));plot_checks<-list()
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
 plot_checks[[j]]<-data.table(Gene=gene,Cells=nrow(q$data),UniqueCells=length(unique(ix)),
   Positive=sum(!is.na(q$data$value)),Zero=sum(is.na(q$data$value)),TitleColor=q$theme$strip.text.x$colour,TitleFace=q$theme$strip.text.x$face)
 plots[[j]]<-ggplotGrob(q);
 if(j==1L) {cairo_pdf(file.path(out,'first_panel_preview.pdf'),width=40/25.4,height=45/25.4,family='Arial');grid::grid.draw(plots[[j]]);dev.off()}
 message(Sys.time(),' rendered scop ',gene)
 rm(q);gc(FALSE)
}
fwrite(rbindlist(plot_checks),file.path(out,'native_plot_checks.tsv'),sep='\t')
shared_heights<-do.call(grid::unit.pmax,lapply(plots,function(g)g$heights));plots<-lapply(plots,function(g){g$heights<-shared_heights;g})
p<-wrap_plots(c(lapply(plots,wrap_elements),rep(list(plot_spacer()),3)),ncol=6) &
  theme(plot.margin=margin(0,0,0,0,'mm'))
g<-patchworkGrob(p)
cairo_pdf(file.path(out,'figS4.pdf'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,family='Arial');grid::grid.draw(g);dev.off()
dev.off();unlink(file.path(out,'measurement.pdf'))
write_json(list(state='COMPLETE',renderer='scop::FeatureDimPlot',expression_palette='Spectral',title_colors=title_colors,feature_headings='italic; colored by marker group',celltype_subtitles=FALSE,cells=2602031,datasets=22,features=features,marker_groups=spec$Marker_Group,layout_columns=6,layout_rows=6,axes=FALSE,raster_pixels=rep(raster_pixels,2),raster_point_radius=raster_point_radius,legend_label_margin_pt=1.5,legend_linear_scale=legend_scale,legend_area_scale=legend_scale^2,figure_width_mm=figure_width_mm,figure_height_mm=figure_height_mm,normalization='log1p(10000 * observed count / full 24659-gene library size)',sampling=FALSE,simulated_expression=FALSE,imputation=FALSE,expression_clipping=FALSE,reclustering=FALSE,annotation_changes=FALSE,each_panel_all_cells=TRUE,positive_cells_drawn_over_zero_cells=TRUE,scale='independent full observed range per gene; see each colorbar',cluster_gene_checks=nrow(check),max_mean_difference=max(abs(check$MeanLog-check$ReferenceMean)),max_detection_difference=max(abs(check$PctPositive-check$ReferencePct)),coordinates_sha256=digest(file=file.path(root,'00_input_audit/compact/embedding_umap.rpca.rds'),algo='sha256')),file.path(out,'figure_audit.json'),auto_unbox=TRUE,pretty=TRUE)
message(Sys.time(),' COMPLETE')
