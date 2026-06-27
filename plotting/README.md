# Integrated atlas plotting

This directory contains the 2.2M-cell integrated human brain atlas plotting
workflow adapted from `HARNexus/code/plotting`.

Run scripts from the repository root:

```sh
bash 04_datasets_plotting.sh
```

Use `bash 04_datasets_plotting.sh T` to force regeneration of existing plots.

Expected inputs:

- `../../data/BrainOmicsData/integration/objects_celltype_plot.rds`
- `../../data/BrainOmicsData/integration/lisi_data.rds`

Outputs:

- `figures/datasets/development_stage_annotation.pdf`
- `figures/datasets/dot_plot_markergenes.pdf`
- `figures/datasets/feature_plots_rpca.pdf`
- `figures/datasets/dim_plots_seurat_clusters_rpca.pdf`
- `figures/datasets/brainregion_stage_rpca.pdf`
- `figures/datasets/dim_plots_datasets_stage_celltype.pdf`
- `figures/datasets/group_heatmap_markergenes.pdf`
- `results/lisi/lisi_results.rds`
- `figures/lisi/lisi_plot_3.pdf`
- `figures/lisi/lisi_plot_2.pdf`
