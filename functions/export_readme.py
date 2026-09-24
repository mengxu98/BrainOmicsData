#!/usr/bin/env python3
"""Write the released package README.

Usage: export_readme.py PACKAGE_DIR
"""
import pathlib, sys

README = '''# A human brain single-cell and single-nucleus transcriptomic resource across age intervals

The package contains 2,602,031 cells or nuclei, 24,659 genes, 22 source datasets, 75 RPCA clusters (C0-C74) and 12 cell types.

Counts are split by source dataset for selective loading and to stay within sparse-matrix index limits. Every shard uses the same ordered gene set.

- `expression/shard_manifest.tsv` indexes the 22 count-matrix shards in 10x Genomics Matrix Market format.
- `metadata/metadata.tsv.gz` contains 21 columns; `metadata/metadata_dictionary.tsv` defines each field.
- `Original_CellType` holds the source-study label; `CellType` holds the adopted 12-type annotation.
- `Age` keeps the reported units and interval precision; `AgeRange` gives the assigned developmental interval. Donor, specimen and library identifiers occupy separate columns; missing libraries are blank.
- `n_counts` and `n_genes` describe the 24,659-gene panel. `percent_mito` is computed from each source count matrix; `metadata/metadata_dictionary.tsv` defines the formula and blank cases.
- `embeddings/` holds 50-dimensional coordinates and matching 2D UMAPs for Raw PCA (`unintegrated_pca`, `unintegrated_umap`), scVI (`scvi`, `scvi_umap`), Harmony (`harmony`, `harmony_umap`) and RPCA (`integrated_pca`, `integrated_umap`).
- `validation/lisi.tsv.gz` holds per-cell dataset iLISI and source-annotation cLISI for the four 50-dimensional representations.
- `provenance/dataset_manifest.tsv` holds one row per source dataset with counts, publication identifiers, processed-input routes, licence terms and the source column behind `Original_CellType`.
- `scripts/readers.zip` contains the R reader (`read_seurat.R`) and the Python reader (`read_h5ad.py`); the Python reader is delivered inside the archive because the deposit endpoint rejects `.py` filenames.
- `provenance/file_manifest.tsv` lists payload dimensions, package schema version, licence, size, SHA256 and MD5. Dimensions are rows x columns (matrix shards: features x cells; embeddings and LISI: cells x coordinates or metrics). Expression-shard licences follow `provenance/dataset_manifest.tsv`; other package files use MIT. `md5sum.txt` also covers the manifest.
'''


def main():
    package = pathlib.Path(sys.argv[1]).resolve()
    (package/'README.md').write_text(README)
    print('export_readme: README.md written')


if __name__ == '__main__':
    main()
