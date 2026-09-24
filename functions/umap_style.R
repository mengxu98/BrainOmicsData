# Shared UMAP style: data panel slightly taller than the full twelve-type legend.
# Based on the four-method panel style; no title, subtitle or audit caption.
umap_style <- list(panel_to_legend_height_ratio=1.12,grid_width_mm=190,grid_height_mm=170,point_size=0.8,raster_dpi=c(1600,1600))
umap_axis_args <- list(lab_size=7,axis_lwd=0.6,xlen_npc=0.10,ylen_npc=0.13,
  arrow_len=0.012,text=element_text(family='Arial'))
umap_theme <- theme(text=element_text(family='Arial',size=7),
  legend.text=element_text(size=6.5),legend.title=element_text(size=7),
  legend.key.height=grid::unit(2.6,'mm'),legend.key.width=grid::unit(2.5,'mm'),
  legend.spacing.x=grid::unit(1,'mm'),legend.box.spacing=grid::unit(2,'mm'),
  plot.margin=margin(3,3,8,4,'mm'))
