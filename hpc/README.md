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
SCIENCEDB_OUT_DIR=/path/to/hpc/home/data/BrainOmicsData/ScienceDB
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

## Harmonized Metadata Schema Maintenance

Run this only when existing RDS metadata columns need to be updated without
rerunning integration:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
sbatch hpc/hpc_metadata.sbatch
```

The maintenance job now updates the full harmonized metadata schema, including
age intervals, raw and standardized sequencing modality, sequencing technology,
reported-age numeric fields and standardized sex labels. Backups are written to
`$INTEGRATION_DIR/schema_backups` by default.

## Clean ScienceDB Matrix Package

The compact reuse package writes directly to
`/path/to/hpc/home/data/BrainOmicsData/ScienceDB` and contains a
10X-compatible matrix, minimal corrected metadata, PCA/UMAP coordinates, reader
scripts and a provenance folder.

Submit after the harmonized metadata schema has been updated:

```sh
cd /path/to/hpc/home/repositories/BrainOmicsData
sbatch hpc/hpc_sciencedb.sbatch
```

The default output layout is:

```text
ScienceDB/
  expression/matrix.mtx.gz
  expression/features.tsv.gz
  expression/barcodes.tsv.gz
  metadata/metadata.tsv.gz
  metadata/dataset_summary.tsv
  metadata/sample_schema_audit.tsv
  embeddings/integrated_pca.tsv.gz
  embeddings/integrated_umap.tsv.gz
  embeddings/unintegrated_umap.tsv.gz
  validation/lisi.tsv.gz
  objects/objects_celltype_plot.rds
  scripts/read_seurat.R
  scripts/read_h5ad.py
  provenance/references.tsv
  README.md
  file_manifest.tsv
  md5sum.txt
```
