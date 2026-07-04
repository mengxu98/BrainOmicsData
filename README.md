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

ScienceDB export now writes the compact reuse package to
`../../data/BrainOmicsData/ScienceDB`. The package contains a 10X-compatible
count matrix, corrected minimal metadata, PCA/UMAP coordinates, reader scripts
and a `provenance/` folder for source/access/citation audit files.

HPC HPC connection details and Slurm entry points are documented in
`hpc/README.md`.
