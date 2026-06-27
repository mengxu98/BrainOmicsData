# HPC HPC jobs

These scripts are the canonical HPC entry points for heavy BrainOmicsData
data preparation. Keep these files synchronized with the local repository before
submitting jobs.

HPC is used for full data integration, batch correction, and ScienceDB data
preparation steps that must read large RDS objects.

## Connection

Web portal:

```text
https://www.scnet.cn/ui/console/index.html#/job-submit
```

SSH:

```sh
ssh -i /Users/mx/Downloads/user_xh5.hpccube.com_RsaKeyExpireTime_2026-09-23_22-40-55.txt \
  -p 22 \
  -o IdentitiesOnly=yes \
  user@xh5.hpccube.com
```

## Paths

Default HPC paths:

```sh
REPO_DIR=/path/to/hpc/home/repositories/BrainOmicsData
INTEGRATION_DIR=/path/to/hpc/home/data/BrainOmicsData/integration
SCIENCEDB_OUT_DIR=/path/to/hpc/home/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource
```

## Data Integration

Submit the full integration job from a HPC terminal:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
sbatch hpc/hpc_integration.sbatch
```

This job runs `03_datasets_integration.sh`, including metadata
harmonization, object merging, PCA/clustering, RPCA/Harmony batch correction,
LISI input generation, and cell-type annotation.

To force all integration steps to rerun:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
OVERWRITE=T sbatch hpc/hpc_integration.sbatch
```

## Age-Interval Schema Maintenance

Run this only when existing RDS metadata columns need to be updated without
rerunning integration:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
sbatch hpc/hpc_add_age_interval_schema.sbatch
```

## ScienceDB Data Preparation

Submit a portal job at the URL above and set the command to run directly:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
bash hpc/hpc_prepare_sciencedb.sh
```

There is no `hpc/hpc_sciencedb_package.sbatch` wrapper script. Do not submit
a nested `sbatch` command from inside a portal job; that creates a short wrapper
job plus a second real job.

For a lightweight metadata-only ScienceDB staging run, use:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
SKIP_HEAVY=1 bash hpc/hpc_prepare_sciencedb.sh
```

Default inputs and output:

```sh
INTEGRATION_DIR=/path/to/hpc/home/data/BrainOmicsData/integration
OBJECT_FILE=/path/to/hpc/home/data/BrainOmicsData/integration/objects_celltypes.rds
SCIENCEDB_OUT_DIR=/path/to/hpc/home/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource
```

## H5AD Conversion

Prepare the conversion environment once:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
bash hpc/hpc_setup_h5ad_env.sh
```

The default setup installs a lightweight `scop` package that provides
`scop::srt_to_h5ad()` without compiling the full `scop` C++ code. To attempt a
full `scop` source installation instead, set `INSTALL_FULL_SCOP=1`.
The conversion script exports Seurat v5 assay layers matching `counts` by
stacking sparse matrices directly, so split layers such as `counts.1`,
`counts.2`, ... do not need to be joined into one large in-memory R matrix.

Test conversion on the smaller marker-plot object first:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
ENV_NAME=brainomics-h5ad \
OBJECT_FILE=/path/to/hpc/home/data/BrainOmicsData/integration/objects_celltype_plot.rds \
OUT_FILE=/path/to/hpc/home/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource/04_processed_objects/human_brain_marker_plot_test.h5ad \
bash hpc/hpc_convert_seurat_to_h5ad.sh
```

Run the full integrated-object conversion in a portal job with sufficient memory:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
ENV_NAME=brainomics-h5ad \
OBJECT_FILE=/path/to/hpc/home/data/BrainOmicsData/integration/objects_celltypes.rds \
OUT_FILE=/path/to/hpc/home/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource/04_processed_objects/human_brain_integrated_full.h5ad \
bash hpc/hpc_convert_seurat_to_h5ad.sh
```

After writing new h5ad files, refresh the ScienceDB manifest:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
/path/to/hpc/home/miniforge3/envs/brainomics-h5ad/bin/Rscript sciencedb/write_sciencedb_manifest.R \
  --out-dir /path/to/hpc/home/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource
```
