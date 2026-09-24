# HPC execution

Local runs and Slurm jobs use the same R 4.5.1 library and Python environment.
Restore it with `environment/restore.sh`, and set
`BRAINOMICS_ENV_ROOT` to that prefix on every node. No site-specific R library,
private Python wheelhouse or cluster module is selected by the analysis scripts.

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
The RPCA job defaults to `1T` memory; set `BRAINOMICS_RPCA_MEMORY` to an appropriate
site value. Other CPU and GPU resource requests are recorded in `hpc/*.sbatch`.
They describe requests, not measured minimum hardware requirements.

Stage 01 supports the Slurm environment installation wrappers. Use local execution
for stages 02–04 and 10–11. Stage 03's `hpc/run_dataset.sh` is a source dispatcher
that runs in place. Stage 12 is an independent upload operation and is excluded
from the default workflow.

`hpc/pull_source_file.sh` copies individual source files over SSH using the host,
port and remote data root in `hpc/local.env`; it is separate from analysis execution.
