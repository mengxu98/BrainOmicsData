# HPC execution

## Site configuration

Copy `hpc/local.env.example` to `hpc/local.env` and set the cluster paths,
partitions and log directory. Set the GPU partition if needed. The numbered
scripts load this file automatically.

Follow the [environment setup](../environment/REPRODUCIBILITY.md#environment).
Set `BRAINOMICS_ENV_ROOT` to a prefix shared by the compute nodes and place the
[required inputs](../environment/REPRODUCIBILITY.md#workflow-dependencies) in the configured paths.

## Submit jobs

Run from the repository root. For integration with prepared inputs:

```sh
BRAINOMICS_EXECUTOR=slurm bash run_pipeline.sh --from 05 --to 09
```

Stage 01 can submit environment installation. Stages 03, 04, 10 and 11 use the
CPU partition. The drivers wait for completion and report job failures. Inside
an existing allocation, run with `BRAINOMICS_EXECUTOR=local`.

To preprocess one dataset and replace its existing outputs:

```sh
BRAINOMICS_EXECUTOR=slurm bash 03_datasets_preprocessing.sh T GSE67835
```

## Job resources and file transfer

For stages 03, 04, 10 and 11, set `BRAINOMICS_CPU_CPUS`, `BRAINOMICS_CPU_TIME`
and `BRAINOMICS_CPU_MEMORY` as needed. Memory falls back to
`BRAINOMICS_SLURM_MEMORY` if `BRAINOMICS_CPU_MEMORY` is empty. If both are
empty, Slurm applies the site default.

Use `BRAINOMICS_RPCA_MEMORY` for RPCA. `BRAINOMICS_SLURM_MEMORY` can also
override the environment-install memory request. See `hpc/*.sbatch` for
other job settings.

`hpc/pull_source_file.sh` transfers individual files over SSH using
`hpc/local.env`.
