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

## Analysis preparation

The numbered workflow and the standalone analysis commands have separate
entry points. Stage 10 consumes prepared metadata, annotations, embeddings and
summary tables. Prepare the required outputs before running that stage; the
numbered workflow does not automatically execute every analysis below.
The figure table in [Figures](#figures) lists the inputs for each plot.

The following sections describe age evaluation, integration sensitivity and
gene-expression reuse. Their commands perform analysis; use the figure entry
points when only redrawing existing results.

### Age evaluation and donor structure

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

### LISI and clustering sensitivity

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

### Gene-expression reuse

This analysis queries the 11 prespecified genes across 22 source datasets and
prepares the expression summaries and paired contrasts used in Figure 5.
It requires the processed source objects, cell-level metadata, adopted cluster
annotation and existing S15 prefrontal-cortex donor pseudobulk inputs.

To extract the genes from every source, aggregate the outputs and draw the
figure, run:

```sh
bash analysis/fig5_rebuild.sh
```

The wrapper runs `analysis/fig5_extract_full_panel.R` once per source, then
`analysis/fig5_prepare_sources.R`, followed by `plotting/fig5.R`. This is a
standalone analysis command, not an automatic part of stage 10. If the four
plotting tables already exist, stage 10 or `Rscript --vanilla plotting/fig5.R`
redraws the figure directly.

Input and output paths:

| Environment variable | Default | Purpose |
|---|---|---|
| `BRAINOMICS_ANALYSIS_DIR` | `results/analysis_run` | Prepared metadata and fixed-context reuse tables |
| `BRAINOMICS_METADATA_FILE` | discovered under the analysis directory | Cell-level metadata |
| `BRAINOMICS_ANNOTATION_TABLE` | `results/annotation/cluster_annotation.tsv` | Adopted cluster annotation |
| `BRAINOMICS_FIG5_FIXED_DIR` | discovered under the analysis directory | S15 prefrontal-cortex donor tables |
| `BRAINOMICS_PROCESSED_DIR` | `<repository>/../../data/BrainOmicsData/processed` | Processed objects for all sources |
| `BRAINOMICS_FIG5_PANEL_DIR` | `results/fig5_full_panel` | Per-dataset 11-gene extraction and source coverage |
| `BRAINOMICS_FIG5_SOURCE_DIR` | `figures` | Summary and four plotting tables |
| `BRAINOMICS_FIG5_OUTPUT_DIR` | `figures` | Panel PDFs and assembled figure |

`analysis/fig5_extract_full_panel.R` reads each source's raw RNA counts and canonical
feature crosswalk. It groups cells by source, canonical donor, age interval,
standardized brain region and adopted cell type, retaining unmeasured genes as
missing rather than zero. Every source cell must be covered exactly once.
Groups require at least 20 cells and 1,000 library counts. The 11 genes are
the prespecified reuse panel in the manuscript.

`analysis/fig5_prepare_sources.R` checks all 22 source-coverage files and 2,602,031 covered
cells. For A it sums pseudobulk within donor and cell type, averages eligible
donors within source, then weights sources equally for each gene/cell type;
each gene is Z-scored across the 12 types. A reports detection separately.
For B it pairs oligodendrocytes and microglia within source, canonical donor,
age interval and brain region. Multiple matched contexts are averaged within
donor, then donors within source. Sources with at least two paired donors
enter the plot. The point/interval summary is the source-equal mean and a
descriptive 95% t interval across source means. The separate two-sided exact
sign test uses the *direction* of each nonzero source mean, with BH correction
across all 11 genes. Its q values are not derived from the t intervals.

For C the prepare script rebuilds paired donor contrasts from the S15
prefrontal-cortex `fixed_panel_donor_counts.tsv.gz` and
`donor_type_eligibility.tsv` inputs. When `paired_ol_micro_donors.tsv` is present, it verifies the
recomputed 440 contrasts and copies the original table formatting. The donor
pseudobulk counts remain a required upstream analysis input. C shows 33 ROSMAP
and 7 SomaMut donors, with raw donor values, per-source interquartile ranges
and per-source means. C is descriptive and has no new significance test.

The four inputs to `plotting/fig5.R` are
`fig5_full_gene_type_source.tsv`, `fig5_full_paired_study_source.tsv`,
`fig5_full_direction_statistics.tsv`, and
`fig5_fixed_paired_donor_source.tsv` in the source directory. The extractor
also writes per-source coverage summaries, and the preparer writes intermediate tables.
Source objects and generated analysis tables stay outside Git; the analysis
and plotting scripts are versioned.

## Figures

Stage 10 runs the figure entry points in manuscript order: Figures 1–5,
then Supplementary Figures S1–S4. It uses prepared analysis outputs; it does
not rerun integration or statistical evaluation. Each entry point exports a
numbered PDF, PNG and LZW-compressed TIFF under `figures/`. Raster resolution
is selected for the manuscript placement width by `functions/export_png.R`.
Only the light and dark Figure 1 SVGs are tracked for GitHub display.

`BRAINOMICS_ANALYSIS_DIR` selects the analysis tree. The other input overrides
are listed in the repository README.

| Figure | Entry point | Panels in reading order | Prepared inputs |
|---|---|---|---|
| 1 | `plotting/fig1.R` | A: resource composition; B: reported age coverage; C: processing workflow | Resource summary and age-coverage tables |
| 2 | `plotting/fig2.R` | A: four-method UMAPs; B: dataset iLISI; C: source-label cLISI; D: donor geometry; E: age prediction; F: RPCA clusters and cell types | Four saved UMAP embeddings, evaluation tables, metadata and annotation |
| 3 | `plotting/fig3.R` | A: cluster markers; B: cell-type markers; C: source-label concordance | Marker evidence, annotation and source-concordance tables |
| 4 | `plotting/fig4.R` | A: original query labels; B: transferred labels; C: label correspondence; D: donor agreement | Prepared query-mapping tables under `tables/egad_mapping/` |
| 5 | `plotting/fig5.R` | A: gene expression and detection; B: paired source contrasts; C: S15 prefrontal-cortex donor contrasts | Four prepared tables from [gene-expression reuse](#gene-expression-reuse) |
| S1 | `plotting/figS1.R` | Age-interval UMAP | Cell metadata, age summaries and RPCA UMAP |
| S2 | `plotting/figS2.R` | Brain-region UMAP | Cell metadata, region summaries and RPCA UMAP |
| S3 | `plotting/figS3.R --redraw` | Separate UMAPs for the 12 cell types | Metadata, annotation and RPCA UMAP |
| S4 | `plotting/figS4.R` | Feature plots for the 33 markers | Retained marker counts, library sizes, marker list, metadata and RPCA UMAP |

### Figure 2: panels A–F

- **A:** `plotting/fig2a.R` draws dataset-labelled UMAPs in Raw, scVI,
  Harmony and RPCA order, each with 2,602,031 cells and a 22-dataset legend.
  It reads `embedding_umap.unintegrated.rds`, `embedding_umap.scvi.rds`,
  `embedding_umap.harmony.rds` and `embedding_umap.rpca.rds`, plus
  `core_metadata_minimal.tsv.gz`. Each coordinate matrix must contain exactly
  the same ordered cell IDs as the metadata, two columns and finite values.
  These are existing analysis outputs; the script does not fit UMAP.
- **B:** dataset iLISI uses the `latent50` / `Dataset` rows of
  `tables/full_lisi_dataset_summary.tsv`.
- **C:** source-label cLISI uses the `latent50` / `Source_Full` rows of the
  same table. Each scheme contains 22 sources × 4 methods. This table is
  produced upstream by `analysis/sensitivity/metrics_finalize_lisi.R`.
- **D:** donor geometry uses `tables/age_signal/donor_structure_by_dataset.tsv`.
- **E:** age prediction uses `tables/age_signal/age_metrics_summary.tsv` and
  `tables/age_signal/age_metrics_by_dataset.tsv`.
- **F:** `plotting/fig2_umap_panels.R` draws the RPCA cluster UMAP on the left
  and cell-type UMAP on the right. Both plots belong to panel F; there is no
  panel G. Cluster legends follow numeric order, C0–C74; cell types follow
  the shared order in `functions/utils.R`.

The statistical tables are resolved under `BRAINOMICS_FIGURE_DATA_DIR` and
checked before any Figure 2 panel is drawn or written. The full LISI summary
is required; an older UMAP summary cannot substitute for it.

The four coordinate files are resolved under `BRAINOMICS_ANALYSIS_DIR`, or
selected explicitly with `BRAINOMICS_RAW_UMAP_FILE`,
`BRAINOMICS_SCVI_UMAP_FILE`, `BRAINOMICS_HARMONY_UMAP_FILE` and
`BRAINOMICS_RPCA_UMAP_FILE`. Stage 10 redraws all panels A–F; a pre-existing
`fig2a.pdf` is no longer required. A fresh clone still needs the prepared
analysis inputs, which are outside version control.

### Supplementary figure redraws

S1 and S2 redraw from metadata and the RPCA embedding. Stage 10 passes
`--redraw` to S3; without that option, its standalone entry point only exports
raster images from the existing PDF. S4 draws all 33 markers from retained
counts, using `functions/figS4_source.R` and `functions/umap_style.R`.

Figure rendering requires the corresponding prepared inputs and plotting
packages. Missing inputs must be supplied before drawing; the figure stage
does not launch upstream analysis to replace them.
