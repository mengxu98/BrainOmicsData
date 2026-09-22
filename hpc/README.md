# HPC helpers

Submission and transfer helpers for the HPC cluster.

```sh
sbatch hpc/hpc_integration.sbatch       # stages 05-06 checkpoint work
sbatch hpc/hpc_parallel_integration.sbatch STAGE=pca RUN_ROOT=<run>
sbatch hpc/hpc_integration_finalize.sbatch
sbatch hpc/hpc_scvi.sbatch              # stage 07
sbatch hpc/hpc_install_scvi.sbatch
sbatch hpc/hpc_install_integration_dependencies.sbatch
```

`run_dataset.sh` dispatches one dataset to the HPC backends used by stage 03.
`pull_source_file_from_hpc.sh` copies a single file from the cluster to the
storage host over SSH.

Stages 01-12 select these jobs through `BRAINOMICS_EXECUTOR=hpc`; the local
executor runs the same scripts in place.
