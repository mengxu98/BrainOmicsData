# Integrated atlas workflow

This directory contains the complete 2.2M-cell human brain atlas integration
workflow adapted from `HARNexus/code/datasets`.

Run scripts from the repository root in this order:

```sh
bash 03_datasets_integration.sh
```

Use `bash 03_datasets_integration.sh T` to force all integration steps to run
again.

Expected processed inputs:

- `../../data/BrainOmicsData/processed/<dataset>/<dataset>_processed.rds`

Major outputs:

- `../../data/BrainOmicsData/integration/objects_list_raw.rds`
- `../../data/BrainOmicsData/integration/metadata_filtered.rds`
- `../../data/BrainOmicsData/integration/objects_list_processed.rds`
- `../../data/BrainOmicsData/integration/objects_raw.rds`
- `../../data/BrainOmicsData/integration/objects_filtered.rds`
- `../../data/BrainOmicsData/integration/objects_integrated.rds`
- `../../data/BrainOmicsData/integration/lisi_data.rds`
- `../../data/BrainOmicsData/integration/objects_celltypes.rds`
- `../../data/BrainOmicsData/integration/objects_celltype_plot.rds`

Workflow summary:

1. `datasets_integration_01.R` loads all processed datasets, harmonizes age
   intervals and brain-region labels, filters incomplete metadata, and saves
   `metadata_filtered.rds`.
2. `datasets_integration_02.R` restricts each dataset to retained cells, merges
   the full object, computes mitochondrial percentage, and removes technical or
   low-information gene categories.
3. `datasets_integration_03.R` runs normalization, HVG selection, PCA,
   clustering, unintegrated UMAP, RPCA integration, Harmony integration, UMAPs,
   and exports `lisi_data.rds`.
4. `datasets_annotation.R` maps 120 Seurat clusters to 9 major cell types and
   exports the full annotated object plus the marker-gene plotting object.
