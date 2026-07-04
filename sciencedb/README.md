# ScienceDB export

This folder exports the ScienceDB package for **An integrated single-cell and
single-nucleus transcriptomic dataset of the human brain across age intervals**.
The export writes directly reusable analysis files rather than the previous
Data-Paper-style archive layout.

Run from the repository root:

```sh
bash 05_sciencedb.sh
```

The default command now performs both steps needed for a public ScienceDB
release: compact package export followed by irreversible public-identifier
anonymization and manifest/MD5 refresh. To refresh anonymization and checksums
for an existing package without regenerating the expression matrix, run:

```sh
bash 05_sciencedb.sh --skip-export
```

Default output:

```text
../../data/BrainOmicsData/ScienceDB
```

Output layout:

```text
ScienceDB/
  expression/matrix.mtx.gz
  expression/features.tsv.gz
  expression/barcodes.tsv.gz
  metadata/metadata.tsv.gz
  metadata/dataset_summary.tsv
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

The expression matrix is exported from Seurat `counts.*` layers as a
10X-compatible Matrix Market package. `metadata/metadata.tsv.gz` stores the
minimal corrected cell metadata used for reuse. `objects/objects_celltype_plot.rds`
stores a lightweight plotting and metadata-inspection Seurat object.
`validation/lisi.tsv.gz` stores per-cell dataset-label LISI values for raw and
RPCA embeddings. `provenance/references.tsv` stores source dataset references,
access links and citation provenance.

The public package anonymization replaces cell, donor, biological-sample and
library identifiers with package-internal anonymous identifiers. The mapping is
not written to disk. Age and age-interval fields are retained as scientific
variables required for age-stratified reuse.

The legacy Data Paper archive workflow is preserved in
`sciencedb/legacy_data_paper_package/`.
