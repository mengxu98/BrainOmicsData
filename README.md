# BrainOmicsData

This repository collects brain single-nucleus/single-cell multi-omics sequencing datasets and provides a workflow-oriented processing pipeline.

Run the 2.2M-cell human brain atlas workflow from the repository root:

```sh
bash 01_datasets_download.sh
bash 02_datasets_preprocessing.sh
bash 03_datasets_integration.sh
bash 04_datasets_plotting.sh
bash 05_sciencedb.sh
```

Pass `T` as the first argument to force rerun steps whose outputs already
exist:

```sh
bash 03_datasets_integration.sh T
bash 04_datasets_plotting.sh T
```

ScienceDB export defaults to the full package, including integrated-object
derived cell metadata, PCA/UMAP coordinates, cluster labels and marker
summaries. Use `SKIP_HEAVY=1 bash 05_sciencedb.sh` for a lightweight
metadata-only staging run.

Shuguang HPC connection details and Slurm entry points are documented in
`hpc/README.md`.
