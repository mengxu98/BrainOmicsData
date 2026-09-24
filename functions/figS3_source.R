#!/usr/bin/env Rscript
# Cell-type locations on the RPCA embedding.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(Matrix);library(Seurat);library(scop);library(patchwork)})
setDTthreads(4);root<-normalizePath(commandArgs(TRUE)[1]);out<-normalizePath(commandArgs(TRUE)[2]);source('functions/utils.R')
a<-fread(file.path(out,'cluster_annotation.tsv'))
if ('CellType' %in% names(a) && !'Working_CellType' %in% names(a)) a[,Working_CellType:=CellType]
classes<-names(brainomics_celltype_colors)
pal<-brainomics_celltype_colors[classes]
a[,Visual_Colour_Class:=Working_CellType]
find_one<-function(filename){candidates<-list.files(root,recursive=TRUE,full.names=TRUE);matches<-candidates[basename(candidates)==filename];if(length(matches)!=1L)stop('Expected one ',filename,' under ',root,'; found ',length(matches));matches[[1L]]}
metadata_file<-Sys.getenv('BRAINOMICS_FIG2_METADATA_FILE',unset='');if(!nzchar(metadata_file))metadata_file<-find_one('core_metadata_minimal.tsv.gz')
embedding_file<-Sys.getenv('BRAINOMICS_RPCA_UMAP_FILE',unset='');if(!nzchar(embedding_file))embedding_file<-find_one('embedding_umap.rpca.rds')
m<-fread(metadata_file,select=c('Cells','Cluster'));e<-readRDS(embedding_file);stopifnot(nrow(m)==2602031L,identical(rownames(e),m$Cells));ii<-match(m$Cluster,a$Cluster);stopifnot(!anyNA(ii),nrow(a)==75L)
pal <- pal[classes]
stopifnot(length(classes)==12L,setequal(classes,a$Working_CellType),all(a$Visual_Colour_Class%in%names(pal)))
fwrite(a[,.(Cells=sum(Cells)),by=Working_CellType],file.path(out,'celltype_counts.tsv'),sep='\t')
md<-data.frame(SpatialGroup=factor(a$Visual_Colour_Class[ii],levels=classes),Cluster=m$Cluster,row.names=m$Cells)
# The empty assay carries plotting metadata; expression is not used here.
carrier<-sparseMatrix(i=integer(),j=integer(),x=numeric(),dims=c(2L,nrow(m)),dimnames=list(c('unusedA','unusedB'),m$Cells));message('Creating plotting container');o<-CreateSeuratObject(counts=carrier,meta.data=md,min.cells=0L,min.features=0L);colnames(e)<-c('UMAP_1','UMAP_2');o[['umap.rpca']]<-CreateDimReducObject(embeddings=e,key='UMAP_',assay='RNA');rm(carrier,md);gc(FALSE)
targets <- c('Differentiating oligodendrocytes','Oligodendrocytes')
focus <- ifelse(o$SpatialGroup %in% targets, as.character(o$SpatialGroup), 'Other cells')
o$OligoLocation <- factor(focus, levels=c(targets,'Other cells'))
focus_pal <- c(pal[targets], 'Other cells'='#DEDEDE')
source('functions/umap_style.R')
base <- umap_theme
make_umap <- function(field,colors,draw_order) {
  p <- scop::CellDimPlot(o,reduction='umap.rpca',group.by=field,
    palcolor=colors,label=FALSE,show_stat=TRUE,seed=11,raster=TRUE,
    raster.dpi=umap_style$raster_dpi,pt.size=umap_style$point_size,xlab='UMAP_1',ylab='UMAP_2',
    legend.title='Cell type',theme_use='theme_blank_axis',
    theme_args=umap_axis_args,combine=FALSE)[[1]]
  # Explicit named breaks/labels preserve heatmap order and the correct counts.
  counts <- table(o[[field,drop=TRUE]])
  cs <- p$scales$get_scales('colour')
  cs$limits <- names(colors);cs$breaks <- names(colors)
  cs$labels <- setNames(paste0(names(colors),'(',as.integer(counts[names(colors)]),')'),names(colors))
  p <- order_dim_plot_cells(p,draw_order)
  stopifnot(nrow(p$data)==ncol(o),identical(rownames(p$data),draw_order))
  p+labs(title=NULL,subtitle=NULL,caption=NULL)+base+
    guides(colour=guide_legend(override.aes=list(size=1.7),ncol=1))
}
# Fix the data-panel size itself, independently of page and legend width.
# The overview panel is 12% taller than its complete twelve-type legend.
mm <- function(x) grid::convertUnit(x,'mm',valueOnly=TRUE)
save_grob <- function(g,name,width,height) {
  cairo_pdf(file.path(out,paste0(name,'.pdf')),width=width/25.4,height=height/25.4,family='Arial')
  grid::grid.newpage();grid::grid.draw(g);dev.off()
  ragg::agg_png(file.path(out,paste0(name,'.png')),width=width,height=height,units='mm',res=300,background='white')
  grid::grid.newpage();grid::grid.draw(g);dev.off()
}
# Build on a Cairo device so legend measurements match the exported font metrics.
cairo_pdf(file.path(out,'layout_measurement.pdf'),width=7,height=5,family='Arial')
p <- make_umap('SpatialGroup',pal,sort(colnames(o)))
g <- ggplotGrob(p)
legend <- g$grobs[[which(g$layout$name=='guide-box-right')]]
legend_height <- mm(sum(legend$heights))
panel <- g$layout[g$layout$name=='panel',]
aspect <- as.numeric(g$heights[panel$t])/as.numeric(g$widths[panel$l])
stopifnot(is.finite(aspect),aspect>0,legend_height>0)
panel_height <- legend_height*umap_style$panel_to_legend_height_ratio
panel_width <- panel_height/aspect
fix_panel <- function(g) {
  z <- g$layout[g$layout$name=='panel',]
  g$heights[z$t] <- grid::unit(panel_height,'mm')
  g$widths[z$l] <- grid::unit(panel_width,'mm')
  g
}
g <- fix_panel(g)
width <- mm(sum(g$widths));height <- mm(sum(g$heights))
message('Overview panel ',round(panel_width,2),' x ',round(panel_height,2),' mm; legend ',round(legend_height,2),' mm; page ',round(width,2),' x ',round(height,2),' mm')
dev.off();unlink(file.path(out,'layout_measurement.pdf'))
save_grob(g,'integration_umap',width,height)
rm(g,p);gc(FALSE)
cairo_pdf(file.path(out,'layout_measurement.pdf'),width=7,height=5,family='Arial')
highlight_order <- colnames(o)[order(focus!='Other cells',colnames(o))]
p <- make_umap('OligoLocation',focus_pal,highlight_order)
g <- fix_panel(ggplotGrob(p))
save_grob(g,'oligodendrocyte_location_umap',max(width,mm(sum(g$widths))),height)
fwrite(data.table(Group=focus)[,.(Cells=.N),by=Group],file.path(out,'oligodendrocyte_location_counts.tsv'),sep='\t')
rm(g,p,o);gc(FALSE)
# Each panel uses the full RPCA coordinates, with the selected cell type drawn last.
labels <- a$Working_CellType[ii]
xlimits <- range(e[,1]);ylimits <- range(e[,2])
plots <- vector('list',length(classes))
counts <- vector('list',length(classes))
for (j in seq_along(classes)) {
  type <- classes[j]; selected <- labels==type
  ord <- c(which(!selected),which(selected))
  d <- data.frame(x=e[ord,1],y=e[ord,2],group=factor(ifelse(selected[ord],type,'Other cells'),levels=c('Other cells',type)))
  stopifnot(nrow(d)==2602031L,sum(selected)==sum(a[Working_CellType==type,Cells]))
  caption <- paste0(paste(strwrap(type,width=29),collapse='\n'),'\n',format(sum(selected),big.mark=',',scientific=FALSE),' cells')
  q <- ggplot(d,aes(x=x,y=y,colour=group))+
    scattermore::geom_scattermore(pointsize=0.8,pixels=c(900,900))+
    scale_colour_manual(values=c('Other cells'='#DEDEDE',pal[type]),guide='none')+
    scale_x_continuous(limits=xlimits,expand=expansion(mult=.04))+
    scale_y_continuous(limits=ylimits,expand=expansion(mult=.04))+
    do.call(theme_blank_axis,modifyList(umap_axis_args,list(lab_size=5.5,axis_lwd=.4)))+
    labs(title=caption,x='UMAP_1',y='UMAP_2')+
    theme(text=element_text(family='Arial'),plot.title=element_text(size=7,hjust=.5,lineheight=1.05),
      plot.margin=margin(3,3,6,6,'mm'),legend.position='none',aspect.ratio=aspect)
  # Convert one full-cohort panel at a time; release its multi-million-row data.
  plots[[j]] <- ggplotGrob(q)
  counts[[j]] <- data.table(CellType=type,Highlighted=sum(selected),Background=sum(!selected),Total=nrow(d))
  message('Built panel ',j,'/12: ',type)
  rm(d,q);gc(FALSE)
}
# Equal title-row heights keep every map aligned despite wrapped type names.
shared_heights <- do.call(grid::unit.pmax,lapply(plots,function(g) g$heights))
plots <- lapply(plots,function(g) {g$heights <- shared_heights;g})
combined <- wrap_plots(lapply(plots,wrap_elements),ncol=4)
grid_grob <- patchworkGrob(combined)
save_grob(grid_grob,'celltype_location_umap',umap_style$grid_width_mm,umap_style$grid_height_mm)
fwrite(rbindlist(counts),file.path(out,'celltype_location_counts.tsv'),sep='\t')
stopifnot(sum(rbindlist(counts)$Highlighted)==2602031L)
dev.off();unlink(file.path(out,'layout_measurement.pdf'))
message('Supplementary Figure S3 written')
