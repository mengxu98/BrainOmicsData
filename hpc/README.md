# HPC execution

Restore the environment with `environment/restore.sh` and set `BRAINOMICS_ENV_ROOT` to
the environment prefix accessible from the compute nodes.

## Site configuration

Copy `hpc/local.env.example` to the ignored `hpc/local.env` and configure the
environment path, data/run paths, CPU and GPU partitions, and log folder.
The numbered drivers, including the separate upload driver, load this file
before selecting the environment.

Submit through the numbered entry points so that their configuration is exported:

```sh
BRAINOMICS_EXECUTOR=slurm bash run_pipeline.sh --from 05 --to 09
```

Stages 05–09 run in dependency order and wait for each job's exit status:
HVG/scVI input, PCA/Harmony, scVI, RPCA, then collection and final evaluation.
Set `BRAINOMICS_RPCA_MEMORY` for the RPCA job according to your cluster and
workload. See `hpc/*.sbatch` for the other job resource settings.

Stage 01 supports the Slurm environment installation wrappers. Stages 03, 04,
10 and 11 also submit their existing drivers to the CPU partition when
`BRAINOMICS_EXECUTOR=slurm`. They wait for completion and propagate job failure;
the submitted body runs in place on the allocated compute node. Configure
`BRAINOMICS_CPU_CPUS`, `BRAINOMICS_CPU_MEMORY` and `BRAINOMICS_CPU_TIME` as
needed. Unset requests use the site defaults. Stage 03 accepts an optional
overwrite flag and dataset name, for example:

```sh
BRAINOMICS_EXECUTOR=slurm bash 03_datasets_preprocessing.sh T GSE67835
```

Place downloaded inputs in the configured data root and run the remaining
stages on the cluster in the locked environment. Stage 12 uploads the package
only when selected explicitly.

`hpc/pull_source_file.sh` copies individual source files over SSH using the host,
port and remote data root in `hpc/local.env`; it is separate from analysis execution.
