# Shuguang HPC jobs

These scripts are the canonical Shuguang entry points for heavy BrainOmicsData
data preparation. Keep these files synchronized with the local repository before
submitting jobs.

Shuguang is used for full data integration, batch correction, and ScienceDB data
preparation steps that must read large RDS objects.

## Connection

Web portal:

```text
https://www.scnet.cn/ui/console/index.html#/job-submit
```

SSH:

```sh
ssh -i /Users/mx/Downloads/mengxu1310_xh5.hpccube.com_RsaKeyExpireTime_2026-09-23_22-40-55.txt \
  -p 65061 \
  -o IdentitiesOnly=yes \
  mengxu1310@xh5.hpccube.com
```

## Paths

Default Shuguang paths:

```sh
REPO_DIR=/work/home/mengxu1310/repositories/BrainOmicsData
INTEGRATION_DIR=/work/home/mengxu1310/data/BrainOmicsData/integration
SCIENCEDB_OUT_DIR=/work/home/mengxu1310/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource
```

## Data Integration

Submit the full integration job from a Shuguang terminal:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
sbatch hpc/shuguang_integration.sbatch
```

This job runs `03_datasets_integration.sh`, including metadata
harmonization, object merging, PCA/clustering, RPCA/Harmony batch correction,
LISI input generation, and cell-type annotation.

To force all integration steps to rerun:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
OVERWRITE=T sbatch hpc/shuguang_integration.sbatch
```

## Age-Interval Schema Maintenance

Run this only when existing RDS metadata columns need to be updated without
rerunning integration:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
sbatch hpc/shuguang_add_age_interval_schema.sbatch
```

## ScienceDB Data Preparation

Submit a portal job at the URL above and set the command to run directly:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
bash hpc/shuguang_prepare_sciencedb.sh
```

There is no `hpc/shuguang_sciencedb_package.sbatch` wrapper script. Do not submit
a nested `sbatch` command from inside a portal job; that creates a short wrapper
job plus a second real job.

For a lightweight metadata-only ScienceDB staging run, use:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
SKIP_HEAVY=1 bash hpc/shuguang_prepare_sciencedb.sh
```

Default inputs and output:

```sh
INTEGRATION_DIR=/work/home/mengxu1310/data/BrainOmicsData/integration
OBJECT_FILE=/work/home/mengxu1310/data/BrainOmicsData/integration/objects_celltypes.rds
SCIENCEDB_OUT_DIR=/work/home/mengxu1310/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource
```

## H5AD Conversion

Prepare the conversion environment once:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
bash hpc/shuguang_setup_h5ad_env.sh
```

The default setup installs a lightweight `scop` package that provides
`scop::srt_to_h5ad()` without compiling the full `scop` C++ code. To attempt a
full `scop` source installation instead, set `INSTALL_FULL_SCOP=1`.
The conversion script exports Seurat v5 assay layers matching `counts` by
stacking sparse matrices directly, so split layers such as `counts.1`,
`counts.2`, ... do not need to be joined into one large in-memory R matrix.

Test conversion on the smaller marker-plot object first:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
ENV_NAME=brainomics-h5ad \
OBJECT_FILE=/work/home/mengxu1310/data/BrainOmicsData/integration/objects_celltype_plot.rds \
OUT_FILE=/work/home/mengxu1310/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource/04_processed_objects/human_brain_marker_plot_test.h5ad \
bash hpc/shuguang_convert_seurat_to_h5ad.sh
```

Run the full integrated-object conversion in a portal job with sufficient memory:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
ENV_NAME=brainomics-h5ad \
OBJECT_FILE=/work/home/mengxu1310/data/BrainOmicsData/integration/objects_celltypes.rds \
OUT_FILE=/work/home/mengxu1310/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource/04_processed_objects/human_brain_integrated_full.h5ad \
bash hpc/shuguang_convert_seurat_to_h5ad.sh
```

After writing new h5ad files, refresh the ScienceDB manifest:

```sh
cd /work/home/mengxu1310/repositories/BrainOmicsData
/work/home/mengxu1310/miniforge3/envs/brainomics-h5ad/bin/Rscript sciencedb/write_sciencedb_manifest.R \
  --out-dir /work/home/mengxu1310/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource
```
