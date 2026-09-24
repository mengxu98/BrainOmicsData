# Reproducing the resource analyses

This repository provides the numbered data-processing workflow, analysis and
figure scripts, environment locks, and retained clustering-sensitivity summaries.
Run commands from the repository root with the required source data and analysis
inputs available. Small-data tests and recalculation from retained summaries
check specific parts of the workflow; a complete atlas reconstruction requires
the upstream processing and integration stages.

The ScienceDB V2 data release and its package schema 2.0.0 are independent of
the code tag. The public `scripts/readers.zip` contains the R and Python data
readers; it is not an archive of all analysis code.

## Environment and small-data verification

The integration dependencies and defaults are recorded in
`environment/r-packages.lock.tsv`, `environment/python-scvi-freeze.txt`, and the
`hpc/` scripts. The local verification profile is
`environment/r-local-validation.lock.tsv`. Use the matching profile; the checker
rejects version mismatches. The integration lock records the environment used
for integration, independently of the environment used for plotting or tests.

```sh
BRAINOMICS_PYTHON=/path/to/python-with-scanpy-anndata-pandas \
  BRAINOMICS_REPRO_OUTPUT=/path/to/new-test-logs \
  bash environment/run_clean_room_fixture.sh integration
```

For the local verification profile, use `local-validation` instead of
`integration`. The wrapper runs `BRAINOMICS_FULL_TESTS=1`, including a synthetic
120-cell, 40-gene, two-source fixture and both reader round trips. The tests check
count values, cell/gene order, metadata, and embeddings. Missing required Python
reader dependencies are an error. Optional scVI import/version verification uses
`BRAINOMICS_VERIFY_SCVI_PYTHON=true` and `BRAINOMICS_SCVI_PYTHON`.

The package-contract check reports payloads missing from a local copy; checking
the contents of the full data package requires those payloads. Successful fixture
tests alone do not establish that the full atlas has been reconstructed.

## Age evaluation and donor structure

`analysis/evaluate_current_age_and_donor_structure.R` reads the four
50-dimensional representations, metadata, and canonical donor identities from
the full analysis tree. The age evaluation uses frontal-cortex excitatory
neurons with exact postnatal ages and at least 20 cells per donor. It leaves one
study's age labels out at a time and excludes shared donors from training.
Donor-distance comparisons use the full donor centroid set.

Run with the analysis root and a scope directory containing
`tables/metadata_formal.rds`:

```sh
Rscript --vanilla analysis/evaluate_current_age_and_donor_structure.R \
  /path/to/analysis_run /path/to/scope
```

To recalculate the predictions and summaries from retained donor centroids, use
a new output directory:

```sh
Rscript --vanilla analysis/replay_age_from_centroids.R \
  /path/to/age_signal/donor_centroids.rds /path/to/new-age-replay
```

This writes predictions, fold coverage, study-level and overall age metrics,
and donor-distance correlations. It uses existing centroids and does not
reconstruct the cell embeddings or centroids.

## LISI and clustering sensitivity

The LISI calculation and finalizer are
`analysis/revision_sources/metrics_full_lisi.R` and
`analysis/revision_sources/metrics_finalize_lisi.R`. The numerical patch is
`hpc/patches/lisi_double_precision.patch`; the pinned source metadata are in
`environment/lisi/LISI_PINNED_SOURCE`. Installation is implemented by
`hpc/install_integration_dependencies.sbatch`. The calculation checks these
source pins, so an unpatched LISI build is not equivalent. The finalizer includes
the Ma-source exclusion and the three-linked-source exclusion analyses.

`analysis/revision_sources/full_graph_stability.R` runs five Louvain clustering
settings, including the baseline, on the same 2,602,031-cell graph.
`analysis/revision_sources/summarize_stability.R` summarizes partition stability
and cell-type agreement after majority matching. The retained baseline has
75 clusters and an adjusted Rand index of 1 against the reference partition;
the four alternatives have 77, 72, 64 and 86 clusters.

The parameter records, overlap counts, per-cluster metrics, and aggregate TSVs
are under `provenance/clustering_sensitivity/`. Its
[README](../provenance/clustering_sensitivity/README.md) gives the five
settings, input requirements, execution commands, and a check of the saved
overlap counts. This sensitivity analysis varies clustering seed and resolution
on a fixed graph; it does not vary integration or graph construction.

`analysis/revision_sources/reconstruct_seed_resolution_sensitivity.R` is a
parameter reconstruction helper that reads the baseline labels directly. Use
`full_graph_stability.R` to reproduce all five clustering runs, including the
baseline.

## Figures

Run figure scripts in a separate checkout/output location when preserving
existing figures. `BRAINOMICS_ANALYSIS_DIR` selects the full retained analysis
tree. The plotting scripts read analysis outputs; they do not recompute the
full atlas.

Figure 5: see [`analysis/fig5_README.md`](../analysis/fig5_README.md).
`BRAINOMICS_FIG5_PANEL_DIR` selects the 22-source extraction used to build the
plotting tables. Redrawing from those tables is separate from extracting the
11-gene panel from the 22 processed source objects.

Figure S3: `functions/figS3_source.R` implements the cell-type drawing workflow
with repository helpers. The original source is
`analysis/revision_sources/plot_celltype_locations.R`.
`Rscript plotting/figS3.R --redraw` draws from the analysis inputs; the default
exports PNG from the retained PDF.

Figure S4: the UMAP style is defined in `functions/umap_style.R`. The plotting
dependency is pinned in `environment/scop-plotting.lock.tsv`, with its source
archive and license in `environment/vendor/`. Use the installed package or set
`SCOP_SOURCE_PATH` explicitly. For example:

```sh
mkdir -p /path/to/temp
tar -xzf environment/vendor/scop-8aec27fc-source.tar.gz -C /path/to/temp
SCOP_SOURCE_PATH=/path/to/temp/scop-8aec27fc \
  Rscript --vanilla plotting/figS4.R
```

Figure rendering requires the corresponding analysis inputs and plotting
packages. Passing the workflow's syntax and fixture tests does not verify a
rendered figure.
