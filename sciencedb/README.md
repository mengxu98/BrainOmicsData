# ScienceDB export workflow

This folder contains the export workflow for the Neuroscience Bulletin Data Paper
package:

`Human brain single-cell and single-nucleus transcriptomic resource across age intervals`

The scripts are designed to run on the data server, where the full
`BrainOmicsData` integration outputs are available.

## Run

The root-level `05_sciencedb.sh` is the canonical entry point:

```bash
bash 05_sciencedb.sh
```

For a lightweight metadata-only staging run:

```bash
SKIP_HEAVY=1 bash 05_sciencedb.sh
```

The heavy step reads `objects_celltypes.rds` and exports final cell type,
cluster, UMAP/PCA and QC tables. The legacy
`sciencedb/run_sciencedb_package.sh` delegates to the root-level entry point.

## Notes

- The package currently treats controlled-access raw human data as source data
  that must be requested from the original repositories. ScienceDB should host
  derived metadata, processed outputs, workflow scripts and permitted annotations.
- Dataset DOI rows marked `needs_source_publication_confirmation` must be resolved
  through Zotero or publisher/database pages before final submission.
