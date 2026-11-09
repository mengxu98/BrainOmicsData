# Reproducing the resource analyses

Run the commands from the repository root. Prepare the source data and the
retained inputs for the steps you plan to run.

## Environment

The target is Linux x86_64 with **R 4.5.1** and **Python 3.12.8**. Package
versions are in `environment/r-packages.lock.tsv` and
`environment/python-scvi-freeze.txt`.

Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).
Add it to `PATH` or set `BRAINOMICS_MICROMAMBA` to its executable. Installation
requires access to conda-forge, CRAN, Bioconductor, GitHub, PyPI and the PyTorch
wheel index. The default environment directory is `.brainomics-env/`.
To use a different directory, run
`export BRAINOMICS_ENV_ROOT=/absolute/path/to/brainomics-env` first.
With either directory, install, activate, and verify the environment:

```sh
bash 01_environment.sh --install
source environment/activate.sh
bash 01_environment.sh --verify
```

For a fresh restore, choose a new environment directory and an empty
`XDG_CACHE_HOME`. The activation script clears inherited library and compiler
search paths so that the project environment supplies the dependencies.
Install Arial and make it available to fontconfig. GPU runs of scVI need a
compatible NVIDIA driver. The environment provides `pigz`, Poppler and TeX
tools. When XeLaTeX or required macros are missing, the installer uses TinyTeX
2023 with the corresponding frozen TeX Live repository.

See the [Slurm guide](../hpc/README.md) for cluster setup and job submission.
If compute nodes need a proxy, set the site's approved `http_proxy` or
`https_proxy` before installing. `PIP_INDEX_URL` can select a Python mirror;
the package locks and PyTorch CUDA index still apply.

## Software checks

After activating the environment, run the default check. Use the second
command to also run the small synthetic tests:

```sh
bash tests/test_workflow.sh
BRAINOMICS_FULL_TESTS=1 bash tests/test_workflow.sh
```

The first command checks scripts, data tables, interfaces and the available
local package. The second adds small synthetic tests of the analysis
algorithms, data export and reader round trips.

## Workflow dependencies

The numbered stages cover environment setup, source processing, matrix assembly,
integration, drawing figures and local export. Run the standalone annotation
and analysis commands below when their outputs are needed. Keep input and
output directories separate when reproducing a step.

### Source data and count matrices

Stage 02 downloads the processed source files, including the Li and Ma inputs.
Stage 03 prepares source objects and metadata under
`$BRAINOMICS_DATA_ROOT/processed/`. Place the data and reference files on the
machine where you will run the workflow.

To download the inputs for Li or Ma separately:

```sh
bash download/Li_et_al_2018.sh
bash download/Ma_et_al_2022.sh
```

The Li script places the two author RData files in `raw/Li_et_al_2018/rawData/`
under `Prenatal_scRNA_seq/` and `Adult_snRNA_seq/`. For Li, stage 03 also reads
the series matrix fetched by `bash download/GSE81475.sh`. The Ma script unpacks
the [BrainSCOPE archive](https://brainscope.gersteinlab.org/integrative_files.html)
to `raw/GSE207334/Ma_Sestan_mat.rds` and reuses an existing RDS. These paths
are relative to `$BRAINOMICS_DATA_ROOT`.

The full stage 04 run also needs five retained files from the reference build.
Prepare them separately before running stage 04, since stages 02 and 03 do not
produce them. Place them at these paths relative to `$BRAINOMICS_RUN_ROOT`:

```text
inputs/plan/selected_gene_panel_34387.tsv
inputs/plan/retained_cells_exploration.tsv.gz
code/source_row_audit_map.tsv.gz
named_gene_panel_policy.json
named_gene_panel_24659.tsv
```

When you run stage 04, it stages selected source bundles that are missing under
`$BRAINOMICS_RUN_ROOT/inputs/processed/`. If those bundles are already in place,
check which required files are present before assembling the matrices:

```sh
bash 04_source_counts_assembly.sh --check-inputs
```

To reuse an existing policy, supply `applied_gene_row_exclusions.tsv`,
`allowed_count_addition.tsv` and `applied_gene_removal_policy.json` together.
Otherwise stage 04 builds per-source rules from observed counts, source
identifiers and the fixed gene references. Observations and rules are written
under `run/source_gene_groups/` and `run/source_gene_policy/`. If you
use a transfer record, supply its matching manifest as well.

Download gene references with `bash 02_datasets_download.sh --references`.
Use `BRAINOMICS_GENE_REFERENCE_DIR` for a custom location. For exact
reproduction, use the original NCBI Gene files because the download URLs
serve updated files. For GSE81475, retain `raw/GSE81475/counts.txt`: stage 04
uses the original `ENSG|symbol` row IDs. It resolves source features with the
fixed references, using stable IDs and the GTF `hgnc_id` field for
GENCODE-to-HGNC links.

Run `bash 04_source_counts_assembly.sh --dataset GSE67835` to rebuild a single
source. Give each source comparison its own `BRAINOMICS_RUN_ROOT`. Set
`BRAINOMICS_REFERENCE_RUN_ROOT` if you want to compare new count blocks with
retained counts on the same named panel.

### Integration and annotation

Run `bash run_pipeline.sh --from 05 --to 09` with the stage 04 matrices.
These stages select highly variable genes, build the common PCA input and
compute the Raw PCA, Harmony, scVI and RPCA representations. Stage 09 collects
the representations, exports annotation inputs and builds the RPCA graph and
partition. Outputs default to `$BRAINOMICS_RUN_ROOT/integration`; use
`BRAINOMICS_RESULTS_DIR` for another directory.

The partition uses 50-dimensional RPCA, 20 Euclidean neighbors from Annoy
with 50 trees, and adaptive Gaussian weights summed with their transpose.
Louvain algorithm 1 uses resolution 2, 10 starts, 10 iterations and seed
20260730. Use a fresh output directory for a new graph and partition. Check a
new partition against the adopted one before reusing the cluster-to-type map.

To expand the adopted annotation, supply the retained partition and the
75-row cluster-to-type table. Write the cell assignments to a new directory:

```sh
export BRAINOMICS_CLUSTER_ASSIGNMENTS=/path/to/retained/cluster_assignments.rds
export BRAINOMICS_ANNOTATION_TABLE=/path/to/adopted/cluster_annotation.tsv
export BRAINOMICS_ANNOTATION_DIR=/path/to/new/annotation
Rscript --vanilla annotation/datasets_annotation.R --build
```

Use the same annotation paths in downstream analyses. The marker summary
script `analysis/cluster_markers.R` uses the common count layers, exported RNA
feature IDs and adopted cell assignments. For downstream analysis, supply the
prepared `metadata_working.rds`, `independent_labels_by_cell.rds` and the
source-identifier sidecar.

## Downstream analyses

Use the commands below for the outputs you need. Each starts from the
corresponding prepared or retained inputs.

### Age evaluation and donor structure

`analysis/evaluate_current_age_and_donor_structure.R` reads four 50-dimensional
representations, metadata and donor IDs. The age evaluation uses frontal-cortex
excitatory neurons with exact postnatal ages and at least 20 cells per donor.
It holds out one study's age labels at a time and excludes shared donors from
training. Donor-distance comparisons use the full centroid set.

Set `BRAINOMICS_METADATA_FILE` to select the metadata, or let the script find
`metadata_working.rds` under the analysis tree. The adopted annotation comes
from `BRAINOMICS_ANNOTATION_TABLE`, which defaults to
`results/annotation/cluster_annotation.tsv`.

```sh
Rscript --vanilla analysis/evaluate_current_age_and_donor_structure.R \
  /path/to/analysis_run /path/to/scope
```

If you already have donor centroids, use the following command to generate
predictions, fold coverage, study-level and overall age metrics, and donor
distance correlations:

```sh
Rscript --vanilla analysis/age_from_centroids.R \
  /path/to/age_signal/donor_centroids.rds /path/to/new-age-summary
```

### LISI and clustering sensitivity

Use the saved embeddings and the matching metadata, source labels and adopted
annotation. Set `BRAINOMICS_ANALYSIS_DIR` to the existing analysis directory
and `BRAINOMICS_FIGURE_DATA_DIR` to the directory for the summary tables.
Create the output directory before running each method in both spaces and
finalizing the tables:

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

Set `BRAINOMICS_SOURCE_LABEL_FILE` if `independent_labels_by_cell.rds` is stored
separately. The finalizer includes the Ma-source exclusion and the
three-linked-source exclusion analyses. The LISI source pin is in
`environment/lisi/LISI_PINNED_SOURCE` and the numerical patch is
`hpc/patches/lisi_double_precision.patch`; both are used by the environment
installer.

`analysis/sensitivity/full_graph_stability.R` runs five Louvain settings on
the same 2,602,031-cell graph. `analysis/sensitivity/summarize_stability.R`
summarizes partition stability and cell-type agreement after majority
matching. See the [clustering-sensitivity records and commands](../provenance/clustering_sensitivity/README.md)
for the settings and the retained overlap counts.

### Restricted query identities

`reference_knn_projection()` in `functions/utils.R` returns coordinates,
`Candidate_Prediction`,
`Support` and `Margin`. To obtain a broad identity label, supply
`query_in_scope` based on independent information about each query cell.
The evaluated classes are neurons, immune cells, vascular cells, fibroblasts
and the oligodendrocyte lineage. `Prediction` is filled for confirmed cells
with a supported candidate. Other cells keep their diagnostic candidates with
a status of `scope_not_confirmed`, `outside_demonstrated_scope` or
`unsupported_reference_assignment`.

For the EGAD evaluation, use the retained `query_predictions.tsv.gz` and the
source-class map. Source annotations determine applicability. Choose a new
output directory:

```sh
Rscript --vanilla integration/evaluate_source_concordance.R \
  /path/to/retained/query_predictions.tsv.gz \
  data/source_class_map.tsv \
  /path/to/new-output
```

The command summarizes assignments by source class and donor, with additional
diagnostic tables by donor and source. The per-cell file
`mapping/restricted_query_predictions.tsv.gz` and the retained query
predictions are stored outside Git and ScienceDB V2.

### Gene-expression reuse

This analysis extracts 11 prespecified genes across 22 sources. It uses
processed source objects, cell metadata, adopted cluster annotations and
retained S15 prefrontal-cortex donor pseudobulk inputs. Run the wrapper to
extract counts and generate the expression summaries and paired contrasts:

```sh
bash analysis/fig5_rebuild.sh
```

| Environment variable | Default | Purpose |
|---|---|---|
| `BRAINOMICS_ANALYSIS_DIR` | `results/analysis_run` | Prepared metadata and reuse tables |
| `BRAINOMICS_METADATA_FILE` | Discovered under the analysis directory | Cell-level metadata |
| `BRAINOMICS_ANNOTATION_TABLE` | `results/annotation/cluster_annotation.tsv` | Adopted cluster annotation |
| `BRAINOMICS_FIG5_FIXED_DIR` | Discovered under the analysis directory | S15 prefrontal-cortex donor tables |
| `BRAINOMICS_PROCESSED_DIR` | `processed/` under the source data root | Processed source objects |
| `BRAINOMICS_FIG5_PANEL_DIR` | `results/fig5_full_panel` | Per-source gene extractions and coverage |
| `BRAINOMICS_FIG5_SOURCE_DIR` | `figures` | Summary tables |
| `BRAINOMICS_FIG5_OUTPUT_DIR` | `figures` | Figure output |

The extraction groups cells by source, donor, age interval, brain region and
adopted cell type. Groups require at least 20 cells and 1,000 library counts.
Unmeasured genes are recorded as missing. For the S15 paired donor analysis,
supply `fixed_panel_donor_counts.tsv.gz` and `donor_type_eligibility.tsv`.
The wrapper writes `fig5_full_gene_type_source.tsv`,
`fig5_full_paired_study_source.tsv`, `fig5_full_direction_statistics.tsv` and
`fig5_fixed_paired_donor_source.tsv`.

## Figures

Run `bash 10_analysis_figures.sh` with the prepared metadata, adopted
annotations, saved embeddings, analysis tables and retained query and marker
results. `BRAINOMICS_ANALYSIS_DIR` selects the analysis tree;
`BRAINOMICS_ANALYSIS_RUN_DIR` and `BRAINOMICS_FIGURE_DATA_DIR` can specify the
metadata and table locations. The figure entry points in [plotting/](../plotting/)
use these inputs directly.

## Local data-package export

Run `bash 11_sciencedb_export.sh` with the expression object, partition,
reconciled metadata, source identifiers and LISI results. Use
`BRAINOMICS_FROZEN_DIR` for retained inputs, `BRAINOMICS_EXPRESSION_OBJECT` to
select a different expression object and `BRAINOMICS_CLUSTER_ANNOTATION` to
select the annotation used in the export. If the LISI files are stored
elsewhere, use `BRAINOMICS_FULL_LISI_DIR` for per-cell results and
`BRAINOMICS_LISI_DIR` for dataset-iLISI checkpoints.

Set `BRAINOMICS_WORK_DIR` and `BRAINOMICS_PACKAGE_DIR` to local working and
output directories. This stage builds a local package. Keep published
ScienceDB files separate from the run.

## Verified execution boundary

The full reference build covered 22 sources, 2,602,031 cells or nuclei and
24,659 genes, from processed source files through integration and
75-cluster partitioning. Targeted verification checked counts for two sources
in the locked environment. Separate tests of the public Li and Ma inputs
reproduced the retained Stage 03 count matrices and key fields in the
standardized metadata. Age summaries and label-transfer evaluation used
retained donor centroids and predictions. Integration results and adopted
annotations were carried forward. The supplementary
ExecutionRecord records the inputs, outputs and scope of each stage.
