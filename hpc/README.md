# HPC helpers

Submission and transfer helpers for the HPC cluster.

```sh
sbatch hpc/integration.sbatch       # stages 05-06 checkpoint work
sbatch hpc/parallel_integration.sbatch STAGE=pca RUN_ROOT=<run>
sbatch hpc/integration_finalize.sbatch
sbatch hpc/scvi.sbatch              # stage 07
sbatch hpc/install_scvi.sbatch
sbatch hpc/install_integration_dependencies.sbatch
```

`run_dataset.sh` dispatches one dataset to the HPC backends used by stage 03.
`pull_source_file.sh` copies a single file from the cluster to the
storage host over SSH.

Stages 01-12 select these jobs through `BRAINOMICS_EXECUTOR=slurm`; the local
executor runs the same scripts in place.

## Site configuration

`hpc/local.env` (git-ignored) holds the account-specific values: `BRAINOMICS_HPC_ROOT`,
`BRAINOMICS_HPC_PARTITION`, `BRAINOMICS_HPC_GPU_PARTITION`, `BRAINOMICS_HPC_LOG_DIR`,
`BRAINOMICS_HPC_HOST`, `BRAINOMICS_HPC_PORT` and `BRAINOMICS_HPC_DATA_ROOT`.
Copy `hpc/local.env.example` and fill it in; the submission helpers and
`pull_source_file.sh` source it automatically.
