# Reproducing the resource analyses

This repository provides the numbered data-processing workflow, analysis and
figure scripts, environment specifications, and clustering-sensitivity
summaries. Run commands from the repository root with the required source data
and analysis inputs available.

The ScienceDB V2 data release and its package schema 2.0.0 are independent of
the code tag. The public `scripts/readers.zip` contains the R and Python data
readers; it is not an archive of all analysis code.

## Verified execution boundary

The full-cohort reference workflow was rerun from
processed source files through clustering for all 22 studies comprising
2,602,031 cells or nuclei. It rebuilt source counts, Raw PCA, scVI, Harmony
and RPCA representations, their UMAP coordinates, and the 75-cluster partition
for the 24,659-gene reference. This run generated new integration results and
assignments.

The locked environment was restored and verified. We reprocessed two source
datasets through count reconstruction. Native-feature and final-panel counts
matched the retained counts for the cells included in the reference. The
integration and adopted annotations were retained.

Five age outputs were recalculated from retained donor centroids. All
1,661,798 retained query predictions were postprocessed to apply the scope
restriction. The original model fits and projections were retained. The
clustering sensitivity summaries use retained overlap counts. Separate software checks exercised small inputs and
the figure and export tools. Figures 2–4 use retained coordinates and results.

The supplementary ExecutionRecord gives the stages, environments, inputs,
outputs and verification scope. Its separate input list records the SHA256
values. The analysis code for the revised manuscript is tagged
[v2.0.1](https://github.com/mengxu98/BrainOmicsData/tree/v2.0.1). The data deposit
is ScienceDB V2 (doi 10.57760/sciencedb.41612).

The retained-input commands below require `query_predictions.tsv.gz` or
`donor_centroids.rds` as specified. These inputs and the generated per-cell
`mapping/restricted_query_predictions.tsv.gz` are held locally outside the
ScienceDB V2 archive and code distribution. The supplement reports the
aggregate restricted results by source class and by donor. The evaluation
command also generates donor-by-source and label-transition details for
diagnostic use.

### Restricted query identities

`reference_knn_projection()` in `functions/utils.R` retains coordinates,
`Candidate_Prediction`, neighbor-vote `Support` and `Margin`. `Prediction`
is withheld by default with status `scope_not_confirmed`. A caller may supply
one logical `query_in_scope` value per cell only when independent information
establishes applicability to the five evaluated broad classes. Do not set this
flag from the same neighbor votes or from the predominance of brain cells in
a mixed query. Missing applicability information remains unconfirmed.

For confirmed applicable cells, supported candidate labels are collapsed to
neurons, immune cells, vascular cells, fibroblasts or the oligodendrocyte
lineage. Candidates without a supported broad counterpart have no identity
and status `unsupported_reference_assignment`. Confirmed out-of-scope cells
have no identity and status `outside_demonstrated_scope`. The original
twelve-class candidates are diagnostic outputs and are not validated
identities. Fine excitatory/inhibitory discrimination is not established by
the broad-neuron comparison.

The full EGAD evaluation uses author source annotations as an explicit
applicability check, separate from candidate prediction. Reproduce it with a
new output directory:

```sh
Rscript --vanilla integration/evaluate_source_concordance.R \
  /path/to/retained/query_predictions.tsv.gz \
  data/source_class_map.tsv \
  /path/to/new-output
```

To draw Figure 4, place the separately retained `query_predictions.tsv.gz`
in the figure input directory documented below. The Figure 4 script derives
the panel C counts and panel D donor concordance from those predictions and
`data/source_class_map.tsv`; no separate aggregate tables are needed.

This retained-prediction calculation is part of the existing source-concordance
entry point. Three file arguments select it; invocation without arguments retains
the original source-concordance and marker-summary calculation. The retained
outputs predate this code consolidation; their original execution hashes
remain historical evidence.

The source-class, donor and donor-by-source tables report unsupported counts,
original candidate errors, errors among retained assignments, and correct
assignments per all comparable cells. They retain all source categories and
the erythrocyte failure. Outside the shared vocabulary, direct broad-class
error rates are unavailable; the known erythrocyte absence is recorded
separately. Non-comparable does not mean biologically absent. Unsupported
rates from this metadata gate are policy outcomes, not learned unknown-cell
detection rates. Performance on unannotated unfamiliar queries remains
unestablished; such queries receive projection and diagnostic candidates,
without identity output.

## Environment

The unified target for future reproduction is Linux x86_64. In this target,
data processing, integration, analysis and plotting use **R 4.5.1**, with package
versions recorded in
`environment/r-packages.lock.tsv`. Python analyses use **Python 3.12.8** and the
versions recorded in `environment/python-scvi-freeze.txt`.

### Restore on a new Linux machine

Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
and ensure it is on `PATH`, or set `BRAINOMICS_MICROMAMBA` to its executable.
An internet connection is required for conda-forge, CRAN, Bioconductor, GitHub,
PyPI, and the official PyTorch wheel repository. Run from the repository root.
The default environment directory is `.brainomics-env/` in the repository.
Set `BRAINOMICS_ENV_ROOT` to install or use the environment elsewhere:

```sh
export BRAINOMICS_ENV_ROOT=/absolute/path/to/brainomics-env
bash 01_environment.sh --install
source environment/activate.sh
bash 01_environment.sh --verify
```

For a clean restore check, use a previously nonexistent environment prefix and
set `XDG_CACHE_HOME` to a separate empty directory so previously compiled R
packages are not reused. Activation removes inherited R/Python library paths
and compiler include/library search paths that can mix host and project
toolchains.

Slurm users can set `BRAINOMICS_EXECUTOR=slurm` and configure the computing
environment in `hpc/local.env`. For the environment-install job,
`BRAINOMICS_SLURM_MEMORY` can override the default memory request to meet the
site's CPU-to-memory policy. If compute nodes lack direct internet access,
configure the site's approved `http_proxy`/`https_proxy` before installation.
`PIP_INDEX_URL` can select a site-accessible Python package mirror; pinned
package versions and the official PyTorch CUDA wheel index remain unchanged.

If downloads run on a separate host, transfer the inputs to the compute
environment before running source preprocessing, count assembly, integration,
annotation, statistics, plotting or export.
`BRAINOMICS_EXECUTOR=slurm` submits stages 03, 04, 10 and 11 to the configured
CPU partition and waits for their completion. Configure the CPU, memory and
time requests through `BRAINOMICS_CPU_CPUS`, `BRAINOMICS_CPU_MEMORY` and
`BRAINOMICS_CPU_TIME`. `BRAINOMICS_EXECUTOR=local` means in-place execution on
the current host, including inside an existing Slurm allocation. Use a new
run and output directory for reproduction.

scVI GPU runs require a compatible NVIDIA driver and an available GPU.
The shared prefix also supplies `pigz`, Poppler and TeX commands used by package
export and figure assembly. If XeLaTeX or its required macros are absent, the
installer provisions TinyTeX 2023 under the environment prefix and uses the
corresponding frozen TeX Live repository. Activation selects those commands
before the conda TeX engines. Arial must be installed on the host and visible to
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

## Workflow dependencies

The numbered drivers do not by themselves rebuild every published result from
downloads. The following dependency order distinguishes executable transitions
from retained scientific inputs. Retained inputs must be supplied explicitly;
their presence does not mean they have been regenerated in the target environment.

| Step | Entry point and generated output | Required input or boundary |
|---|---|---|
| Environment | `01_environment.sh --install`, then `01_environment.sh --verify` | Linux x86_64, micromamba, network access, the documented fonts and system requirements; a compatible NVIDIA driver for GPU training |
| Source processing | `03_datasets_preprocessing.sh`: per-source objects and canonical metadata under `$BRAINOMICS_DATA_ROOT/processed/` | Downloaded processed source data; existing outputs are retained unless overwrite is requested |
| Assembly preparation | Stage 04 reads `$BRAINOMICS_RUN_ROOT/inputs/processed/` | Actual execution stages missing selected stage-03 bundles. Existing snapshots are retained. The fixed cell list, target gene panels and source-row mapping remain explicit scientific inputs |
| Common matrices | `04_source_counts_assembly.sh`: source-gene observations, row policies, cleaned chunks, named-panel source objects and common matrices | Missing row-exclusion and count-addition tables are generated by the count driver from source-count observations, source feature identifiers and fixed gene references. The existing cohort and feature-to-gene assignments are retained |
| Integration | Stages 05–09: HVGs, common PCA, Raw/Harmony/scVI/RPCA representations and latent-space dataset iLISI | Stage-04 matrices; scVI needs the configured accelerator and RPCA needs the documented memory. Stage 09 also exports annotation latent representations and RNA feature identifiers |
| Partition | Stage 09 calls `functions/run_parallel_integration.R` in `cluster` mode to generate `cluster_assignments.rds` and the weighted RPCA graph | Uses the completed RPCA checkpoint and the recovered historical graph/clustering algorithm. A new partition must be checked against the adopted partition before reusing the existing cluster-to-type mapping; the sensitivity script consumes this graph and partition |
| Adopted annotation | `annotation/datasets_annotation.R --build`: cell assignments and type counts | The retained partition plus the adopted 75-row cluster-to-type table. This expands an existing decision; it neither reclusters nor reannotates cells |
| Marker summaries | `analysis/cluster_markers.R`: per-source checkpoints and marker summaries under the integration annotation tree | Common count layers, exported RNA feature identifiers and adopted cell assignments. These outputs do not automatically replace the separately retained evidence tables used by Figure 3 and S4 |
| Analysis metadata and source labels | Retained `metadata_working.rds`, `independent_labels_by_cell.rds` and identifier sidecar | Their complete historical preparation is not connected to the numbered stages. Do not substitute input metadata for the final donor/specimen reconciliation or source-label curation |
| Statistics | Age, full LISI, stability and Figure 5 commands below | The matching metadata, embeddings, source labels, adopted annotation, and any branch-specific retained inputs; these commands are separate from stage 10 |
| Query mapping | Retained `query_predictions.tsv.gz`; `integration/evaluate_source_concordance.R` in three-argument mode generates restricted reports | The retained predictions are required separately. Recomputing these reports does not repeat query preparation, fitting or projection; the zero-argument mode uses an older, different query cohort |
| Figures | `10_analysis_figures.sh` | Prepared statistics, metadata, eight saved embeddings, Figure 1 summaries, Figure 3 marker/concordance evidence, query results and supplementary marker counts. Missing historical summary producers remain retained-input boundaries |
| Local export | `11_sciencedb_export.sh` | Expression object, partition, reconciled metadata, source-identifier sidecar and full LISI outputs. The historical default is the frozen clustered object; `BRAINOMICS_EXPRESSION_OBJECT` can select a newly generated integration object. Export does not create these inputs or rebuild the atlas |

Stage 04 requires the following retained scientific inputs in addition to the
processed source bundles:

```text
inputs/plan/selected_gene_panel_34387.tsv
inputs/plan/retained_cells_exploration.tsv.gz
code/source_row_audit_map.tsv.gz
named_gene_panel_policy.json
named_gene_panel_24659.tsv
```

Actual execution stages missing selected processed bundles under
`inputs/processed/` and retains existing snapshots. The availability-only check
below expects these bundles to be present and does not stage inputs.

Check these paths before starting expensive matrix construction:

```sh
bash 04_source_counts_assembly.sh --check-inputs
```

This checks input-file availability, not source-data identity or successful
execution. The count builder does not require an input-manifest hash or an old
transfer-completion marker. For a local matrix build,
the common-matrix assembler validates each source object directly. If a
transfer record or manifest is supplied, the matching pair remains required
and is checked. These input checks do not relax cell, feature, count or ordering
checks or establish that a scientific stage has completed.

The legacy `applied_gene_row_exclusions.tsv`, `allowed_count_addition.tsv` and
`applied_gene_removal_policy.json` can be supplied together to reuse an existing
policy. When they are absent, the count driver generates source-gene observations
and per-source policies under `run/source_gene_groups/` and
`run/source_gene_policy/`. No final matrix is used to infer the policy. Source
feature IDs and study-reported reference annotations are checked against the
fixed gene-reference inputs. Unresolved identities and unsupported positive
count additions remain excluded.

For GSE81475, keep `raw/GSE81475/counts.txt` with the processed bundle. The
policy builder reads its original `ENSG|symbol` row identifiers, which the
processed object's symbol aggregation does not retain. Those source IDs are
checked against the fixed references; conflicting or unresolved assignments
remain excluded. Historical IDs require support from official reference names;
an alias alone does not establish identity. GENCODE-to-HGNC links use explicit
stable-ID cross-references and the GTF `hgnc_id` field, rather than a shared gene
symbol. This restores the original source-evidence route without changing the
preprocessing counts or the frozen feature-to-gene assignments.

Download the reference inputs with `bash 02_datasets_download.sh --references`.
If downloads run on a separate host, transfer `raw/gene_references/` to the
compute environment. The reference manifest lists file names and source URLs. The NCBI Gene URLs
are updated over time; use the retained NCBI Gene files from the original run
for exact reproduction. Set `BRAINOMICS_GENE_REFERENCE_DIR` if stored elsewhere.

For a source-by-source rebuild, use
`bash 04_source_counts_assembly.sh --dataset GSE67835`. This does not collect
the full common matrix. Use a separate `BRAINOMICS_RUN_ROOT` for each source
comparison so its inputs and generated policies remain separate. To compare
rebuilt counts with retained named-panel source objects, optionally set
`BRAINOMICS_REFERENCE_RUN_ROOT`. When set, every newly cleaned count block is
compared with the corresponding retained counts. A difference stops the source
build before further assembly or integration. Neither generating new policies nor passing
this comparison constitutes a new RPCA or scVI execution.

To expand the adopted partition in an independent annotation output directory:

```sh
export BRAINOMICS_CLUSTER_ASSIGNMENTS=/path/to/retained/cluster_assignments.rds
export BRAINOMICS_ANNOTATION_TABLE=/path/to/adopted/cluster_annotation.tsv
export BRAINOMICS_ANNOTATION_DIR=/path/to/new/annotation
Rscript --vanilla annotation/datasets_annotation.R --build
```

The command refuses to replace an existing assignment file and checks the
2,602,031 cells, 75 cluster sizes and 12 types before writing. Downstream
assignment readers use the same annotation-directory and annotation-table
variables. A newly fitted integration must not be given adopted cluster IDs
solely because it contains the same cells: the partition must be independently
established or explicitly retained and described as such.

### Execution needed before broader reproduction claims

First restore the complete environment into a new prefix and execute the
actual source-processing entry points in an independent data tree. Then recover
and check the assembly policies and stage-03-to-stage-04 input contract. Only
after those dependencies are available can a clean matrix/integration build
be attempted. None of these steps has been demonstrated by the synthetic
fixture or the retained-prediction replay.

The recovered partition algorithm uses 50-dimensional RPCA, Annoy Euclidean
20-neighbor search with 50 trees, adaptive Gaussian directed weights summed
with their transpose, and Louvain algorithm 1 at resolution 2, with 10 starts,
10 iterations and seed 20260730. Stage 09 now connects this calculation to the
existing RPCA output and refuses to overwrite graph or partition outputs.
This connection has not yet been verified by a new full-cohort clustering run.
Its output is not automatically assigned the adopted biological labels.

With retained scientific inputs, real-data downstream replay can instead
start at annotation expansion, age/donor statistics, the eight method/space
LISI calculations and finalization, Figure 5 extraction/aggregation, and figure
rendering. Such a run must identify those retained inputs and must not be
reported as a full build from downloads. The existing guide commands describe
each branch; the compact 120-cell fixture remains a software regression test.

Use separate input, run, work, figure and package paths for any replay. Do not
point an exporter at the published ScienceDB directory or select upload stage
12. Local reproduction does not require changing ScienceDB V2.

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

Run every method in both evaluated spaces before finalizing the tables:

```sh
for space in latent50 umap2; do
  for method in Raw RPCA Harmony scVI; do
    "$BRAINOMICS_RSCRIPT" --vanilla analysis/sensitivity/metrics_full_lisi.R \
      "$BRAINOMICS_ANALYSIS_DIR" "$BRAINOMICS_FIGURE_DATA_DIR" "$method" "$space"
  done
done
"$BRAINOMICS_RSCRIPT" --vanilla analysis/sensitivity/metrics_finalize_lisi.R \
  "$BRAINOMICS_FIGURE_DATA_DIR"
```

Both directories must already exist. This is a full-cohort statistical
calculation from supplied embeddings, not a rerun of integration. Stage 09's
dataset-iLISI checkpoints alone do not produce all of these source-label and
formal-type summaries.

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
| `BRAINOMICS_PROCESSED_DIR` | `processed/` under the source data root | Processed objects for all sources |
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

To reuse the retained Figure 4A/B panels, set `BRAINOMICS_FIG4_REUSE_AB_DIR`
to a directory containing `fig4a.pdf` and `fig4b.pdf` when running
`plotting/fig4.R`. Without this setting, A/B are redrawn from the retained
query predictions. The retained query predictions and source-class map are required in either
mode.

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
