# Reproducing the resource analyses

This repository provides the numbered data-processing workflow, analysis and
figure scripts, environment specifications, and clustering-sensitivity
summaries. Run commands from the repository root with the required source data
and analysis inputs available.

The ScienceDB V2 data release and its package schema 2.0.0 are independent of
the code tag. The public `scripts/readers.zip` contains the R and Python data
readers; it is not an archive of all analysis code.

## Environment

The supported platform is Linux x86_64. Data processing, integration, analysis
and plotting use **R 4.5.1**, with package versions recorded in
`environment/r-packages.lock.tsv`. Python analyses use **Python 3.12.8** and the
versions recorded in `environment/python-scvi-freeze.txt`.

### Restore on a new Linux machine

Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
and ensure it is on `PATH`, or set `BRAINOMICS_MICROMAMBA` to its executable.
An internet connection is required for conda-forge, CRAN, Bioconductor, GitHub,
PyPI, and the official PyTorch wheel repository. Run from the repository root:

```sh
export BRAINOMICS_ENV_ROOT=/absolute/path/to/brainomics-env
bash 01_environment.sh --install
source environment/activate.sh
bash 01_environment.sh --verify
```

Slurm users can set `BRAINOMICS_EXECUTOR=slurm` and configure the computing
environment in `hpc/local.env`.

scVI GPU runs require a compatible NVIDIA driver and an available GPU.
The shared prefix also supplies `pigz`, Poppler and TeX commands used by package
export and figure assembly. Arial must be installed on the host and visible to
fontconfig; the environment check rejects a substitute font.

## Software checks

After restoring and activating the environment, run from the repository root:

```sh
bash tests/test_workflow.sh
BRAINOMICS_FULL_TESTS=1 bash tests/test_workflow.sh
```

The first command checks syntax, agreement between data tables and scripts,
script interfaces, and the available local package contract. The second command
adds small synthetic tests:

- `tests/test_synthetic_pipeline.R` checks age mapping and the correspondence
  between integration methods and LISI outputs. It also exports a dataset with
  120 cells, 40 genes and 2 expression shards and reads it back through the R
  and Python readers. `tests/test_reader_roundtrip.py` compares the Python
  output with the exported counts, cell/gene order, metadata and embeddings.
- `tests/test_integration_algorithms.R` exercises PCA, RPCA, Harmony, neighbor
  construction and clustering on a small synthetic dataset.
- `tests/test_lisi_neighbors.R` checks the LISI neighbor calculations; additional
  tests cover age-interval decisions, query normalization and metadata handling.

Installing the environment or passing `01_environment.sh --verify` checks the
software setup. The commands above execute the software tests. The synthetic
tests do not reconstruct the full atlas or train scVI; full-data analyses
require the source datasets and prepared inputs described in this guide.

## Age evaluation and donor structure

`analysis/evaluate_current_age_and_donor_structure.R` reads the four
50-dimensional representations, metadata, and canonical donor identities from
the full analysis tree. The age evaluation uses frontal-cortex excitatory
neurons with exact postnatal ages and at least 20 cells per donor. It leaves one
study's age labels out at a time and excludes shared donors from training.
Donor-distance comparisons use the full donor centroid set.

Run with the analysis root and an output scope. Set
`BRAINOMICS_METADATA_FILE` to select the metadata explicitly; otherwise the
script locates `metadata_working.rds` under the analysis tree. Cell types are
mapped from the adopted cluster annotation selected by
`BRAINOMICS_ANNOTATION_TABLE` (default:
`results/annotation/cluster_annotation.tsv`):

```sh
Rscript --vanilla analysis/evaluate_current_age_and_donor_structure.R \
  /path/to/analysis_run /path/to/scope
```

To calculate the predictions and summaries from donor centroids, use
a new output directory:

```sh
Rscript --vanilla analysis/age_from_centroids.R \
  /path/to/age_signal/donor_centroids.rds /path/to/new-age-summary
```

This writes predictions, fold coverage, study-level and overall age metrics,
and donor-distance correlations from the supplied centroids.

## LISI and clustering sensitivity

The LISI calculation and finalizer are
`analysis/sensitivity/metrics_full_lisi.R` and
`analysis/sensitivity/metrics_finalize_lisi.R`. The numerical patch is
`hpc/patches/lisi_double_precision.patch`; the pinned source metadata are in
`environment/lisi/LISI_PINNED_SOURCE`. Installation is implemented by
`environment/restore_r.R` as part of the unified environment restore. The
calculation uses these source pins. The finalizer includes
the Ma-source exclusion and the three-linked-source exclusion analyses.

`analysis/sensitivity/full_graph_stability.R` runs five Louvain clustering
settings, including the baseline, on the same 2,602,031-cell graph.
`analysis/sensitivity/summarize_stability.R` summarizes partition stability
and cell-type agreement after majority matching. The baseline has
75 clusters and an adjusted Rand index of 1 against the reference partition;
the four alternatives have 77, 72, 64 and 86 clusters.

The parameter records, overlap counts, per-cluster metrics, and aggregate TSVs
are under `provenance/clustering_sensitivity/`. Its
[README](../provenance/clustering_sensitivity/README.md) gives the five
settings, input requirements, execution commands, and a check of the saved
overlap counts. This sensitivity analysis varies clustering seed and resolution
on a fixed graph; it does not vary integration or graph construction.

## Figures

`BRAINOMICS_ANALYSIS_DIR` selects the analysis tree used by the figure scripts.

Figure 5: see [`analysis/fig5_README.md`](../analysis/fig5_README.md).
`BRAINOMICS_FIG5_PANEL_DIR` selects the 22-source extraction used to build the
plotting tables. Redrawing from those tables is separate from extracting the
11-gene panel from the 22 processed source objects.

Figure S3: `functions/figS3_source.R` implements the cell-type drawing workflow
with repository helpers.
`Rscript plotting/figS3.R --redraw` draws from the analysis inputs; the default
exports PNG from the corresponding PDF.

Figure S4: the UMAP style is defined in `functions/umap_style.R`. Run:

```sh
Rscript --vanilla plotting/figS4.R
```

Figure rendering requires the corresponding analysis inputs and plotting
packages.
