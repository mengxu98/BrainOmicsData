# A human brain single-cell and single-nucleus transcriptomic resource across age intervals

BrainOmicsData is the code repository for this resource.

**2,602,031 cells or nuclei · 24,659 genes · 22 source datasets · 75 clusters · 12 cell types**

Version 2 is published at ScienceDB under the existing DOI [10.57760/sciencedb.41612](https://doi.org/10.57760/sciencedb.41612). The analysis code for the revised manuscript is tagged [v2.0.1](https://github.com/mengxu98/BrainOmicsData/tree/v2.0.1). The ScienceDB files remain unchanged.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/fig1-night.svg">
  <source media="(prefers-color-scheme: light)" srcset="figures/fig1.svg">
  <img alt="Overview of the integrated human brain atlas" src="figures/fig1.svg" width="100%">
</picture>

## Pipeline

The numbered drivers run from the repository root. The default analysis workflow ends at stage 11; stage 12 is a separate upload command.

| Stage | Driver | What it does |
|---|---|---|
| 01 | `01_environment.sh` | Install or verify the analysis environment |
| 02 | `02_datasets_download.sh` | Download each source dataset |
| 03 | `03_datasets_preprocessing.sh` | Reconstruct and standardize each source |
| 04 | `04_source_counts_assembly.sh` | Build cleaned count chunks, assemble the named gene panel, and accept the common matrices |
| 05 | `05_hvg_scvi_input.sh` | Select 3,000 HVGs and export the scVI input |
| 06 | `06_pca_harmony.sh` | Raw PCA, raw UMAP, Harmony |
| 07 | `07_scvi.sh` | scVI latent space |
| 08 | `08_rpca.sh` | RPCA latent space |
| 09 | `09_final_assembly.sh` | Collect completed reductions, calculate the RPCA graph/partition, import scVI, evaluate latent spaces and export annotation representations |
| 10 | `10_analysis_figures.sh` | Main and supplementary figures from the prepared analysis outputs |
| 11 | `11_sciencedb_export.sh` | Build the cell-level QC and annotation inputs when absent, then export and verify a local data package |
| 12 | `12_sciencedb_upload.sh` | Transfer the sealed package |

`bash run_pipeline.sh` runs stages 01–11; `--list`, `--from`, `--to`, and `--only` select subsets. Stage 09 follows completion of all integration methods. Upload is selected explicitly with `--only 12` and is not required for reproduction.

`BRAINOMICS_EXECUTOR=local` (default) runs stages in place. Stages 01 and 03–11 also support `BRAINOMICS_EXECUTOR=slurm`. Data preparation, figure assembly and local package export submit the existing drivers to a CPU allocation; integration uses its method-specific job scripts. All submissions wait for the job's exit status. See [`hpc/README.md`](hpc/README.md).

For HPC reproduction, place the downloaded source data and the retained inputs
listed in the guide on the cluster. Restore and verify the environment with
stage 01 before running `--from 03 --to 11`. Inside a CPU allocation, use
`BRAINOMICS_EXECUTOR=local` for stages without a Slurm wrapper. A source-policy
comparison uses stages 03 and 04; run stages 05–09 for a new integration.

The numbered stages use the downloaded source datasets together with prepared
analysis inputs. Stage 04 uses fixed cell/gene mappings and can regenerate
source-row exclusions and count-addition rules from the source counts and frozen
gene references;
stage 10 reads the analysis metadata, statistical tables and saved UMAP coordinates;
stage 11 reads the expression and annotation trees. Their locations can be set
with the environment variables listed below. Figure 2 uses the saved UMAP
coordinates for all four methods; the panel-to-input mapping is described in
the [reproducibility guide](environment/REPRODUCIBILITY.md#figures).

### Environment

See the [environment setup instructions](environment/REPRODUCIBILITY.md#environment)
for software requirements and installation commands.

## Layout

| Path | Contents |
|---|---|
| `functions/` | Shared helpers, figure assembly and support scripts, package exporters, normalization scripts |
| `integration/` | Integration and evaluation modules |
| `processing/`, `download/` | Per-dataset reconstruction and download scripts |
| `annotation/` | Annotation inputs and adopted-label export |
| `analysis/`, `plotting/` | Analysis summaries and figure-numbered plotting scripts |
| `sciencedb/` | Package metadata, manifest, reader and anonymization utilities |
| `results/` | Local analysis inputs used by the figure scripts (outside version control) |
| `environment/` | Package locks and restore/verification scripts |
| `hpc/` | HPC submission scripts (`hpc/*.sbatch`, site settings in `hpc/local.env`) and transfer helpers |
| `data/` | Source access table, donor crosswalks, feature metadata |
| `tests/` | Workflow tests and package contract checks |

Analysis preparation, input requirements and plotting commands are documented
in the [reproducibility guide](environment/REPRODUCIBILITY.md).
Its [verified execution boundary](environment/REPRODUCIBILITY.md#verified-execution-boundary)
distinguishes retained results, actual recalculations and synthetic checks from
a full clean rebuild. The query restriction withholds identity by
default and only returns evaluated broad classes with independently confirmed
applicability. It does not implement learned unknown-cell detection. See
[restricted query identities](environment/REPRODUCIBILITY.md#restricted-query-identities).
The source-class crosswalk is provided in `data/`. Figure 4C/D summaries are
calculated from the retained query predictions. The predictions and donor
centroids are required separately as described in the reproducibility guide. The submitted supplementary ExecutionRecord worksheet
summarises the actual execution boundary.
[Analysis preparation](environment/REPRODUCIBILITY.md#analysis-preparation)
describes the standalone analyses that supply stage 10;
[Figures](environment/REPRODUCIBILITY.md#figures) lists the drawing entry points
in order, from Figures 1–5 to S1–S4. Each exports PDF, PNG and LZW-compressed
TIFF under `figures/`. Shared plotting helpers live under `functions/`.

The [workflow dependency table](environment/REPRODUCIBILITY.md#workflow-dependencies)
identifies retained inputs and transitions that are not generated by the
numbered stages. In particular, stage 03 does not prepare the frozen stage-04
cell/gene selection inputs, and integration does not recreate the biological annotation,
final metadata reconciliation or every figure summary. Check stage-04 input
availability with `bash 04_source_counts_assembly.sh --check-inputs` before a
matrix rebuild. Fetch gene-reference files with `bash 02_datasets_download.sh --references`.
For exact reproduction, retain the NCBI Gene files used in the original run;
the download URLs provide the current releases. The count driver scans source gene groups and generates its row-policy
tables when the legacy tables are absent. `--dataset NAME` rebuilds a single
source without collecting the full common matrix. To compare regenerated counts
with retained final-panel source objects, optionally set
`BRAINOMICS_REFERENCE_RUN_ROOT`; any count difference stops execution before
integration. The annotation entry point can expand an explicitly supplied
adopted partition with `--build`; this does not infer new cell identities.

Local analysis inputs live under `results/`: `analysis_run/` (figure inputs), `frozen_run/` (stage-11 inputs), `annotation/`, `run_root/` (stage-04 inputs and outputs) and `work/` (stage-11 work directory). `results/`, `figures/` (except the light and dark Figure 1 SVGs) and the manuscript folders are outside version control. PNG and TIFF exports remain local for manuscript preparation; GitHub displays the corresponding SVG for the reader's color theme.

Analysis modules that use the logical name `integration_25` resolve it through
`BRAINOMICS_RESULTS_DIR`. Figures and annotation use the 2,602,031-cell inputs
under `results/analysis_run/` and `results/annotation/`.

## Environment variables

| Variable | Purpose | Used by |
|---|---|---|
| `BRAINOMICS_DATA_ROOT` | Directory for source downloads and processed objects | 02, 03, 11 |
| `BRAINOMICS_RESULTS_DIR` | Integration outputs; defaults to `$BRAINOMICS_RUN_ROOT/integration` | 04–09, 11 |
| `BRAINOMICS_RUN_ROOT` | Input policies and run outputs; defaults to `<repository>/results/run_root` | 04–09, 11 |
| `BRAINOMICS_FROZEN_DIR` | Frozen results used for local package assembly | 11 |
| `BRAINOMICS_EXPRESSION_OBJECT` | Explicit expression/reduction object for local export; defaults to the frozen clustered object | 11 |
| `BRAINOMICS_ANALYSIS_DIR` | Prepared analysis and figure inputs | 10, 11 |
| `BRAINOMICS_ANALYSIS_RUN_DIR` | Directory containing `01_metadata/metadata_working.rds` | 10 |
| `BRAINOMICS_FIGURE_DATA_DIR` | Directory containing the figure summary tables | 10 |
| `BRAINOMICS_REFERENCE_SUMMARY` | Figure 1 resource-summary directory | 10 |
| `BRAINOMICS_FIG2_METADATA_FILE` | Cell, dataset and cluster table used by Figures 2A and 2F | 10 |
| `BRAINOMICS_RAW_UMAP_FILE` | Saved Raw PCA UMAP coordinates for Figure 2A | 10 |
| `BRAINOMICS_SCVI_UMAP_FILE` | Saved scVI UMAP coordinates for Figure 2A | 10 |
| `BRAINOMICS_HARMONY_UMAP_FILE` | Saved Harmony UMAP coordinates for Figure 2A | 10 |
| `BRAINOMICS_RPCA_UMAP_FILE` | RPCA UMAP coordinates used by Figures 2 and S1–S3 | 10 |
| `BRAINOMICS_METADATA_FILE` | Cell-level analysis metadata | 11 and standalone analyses |
| `BRAINOMICS_ANNOTATION_TABLE` | Adopted 75-cluster annotation used by assignment expansion, analyses and figures | annotation, 10 |
| `BRAINOMICS_ANNOTATION_DIR` | Adopted cell assignments and type counts; defaults to `results/annotation` | annotation expansion and assignment readers |
| `BRAINOMICS_CLUSTER_ANNOTATION` | Versioned cluster annotation used for the cell-level package table | 11 |
| `BRAINOMICS_SOURCE_LABEL_FILE` | Harmonized source-label table used for cLISI | standalone LISI analysis |
| `BRAINOMICS_MARKER_COUNTS_DIR` | Full-cell marker-count inputs for Figure S4 | 10 |
| `BRAINOMICS_MARKER_ORDER_FILE` | Ordered marker list for Figure S4 | 10 |
| `BRAINOMICS_MARKER_EVIDENCE_FILE` | Cluster-level expression reference for Figure S4 | 10 |
| `BRAINOMICS_FULL_LISI_DIR` | Per-cell cLISI result directory used during package assembly | 11 |
| `BRAINOMICS_CLUSTER_ASSIGNMENTS` | Adopted partition used for annotation expansion and cell-level QC | annotation, 11 |
| `BRAINOMICS_FIG5_PANEL_DIR` | Per-source Figure 5 extraction tables | standalone Figure 5 rebuild |
| `BRAINOMICS_FIG5_FIXED_DIR` | S15 prefrontal-cortex donor tables for Figure 5C | standalone Figure 5 rebuild |
| `BRAINOMICS_FIG5_SOURCE_DIR` | Figure 5 plotting tables | 10 |
| `BRAINOMICS_FIG5_OUTPUT_DIR` | Figure 5 panel and assembled figures | 10 |
| `BRAINOMICS_PACKAGE_DIR` | Local package output; defaults to `$BRAINOMICS_RUN_ROOT/package` | 11 |
| `BRAINOMICS_WORK_DIR` | Export work directory (cell-level QC and annotation tables) | 11 |
| `BRAINOMICS_LISI_DIR` | Per-cell iLISI checkpoints | 11 |

Stage 11 writes a local package. The published ScienceDB files are not changed by stages 01–11.

## Checks

```sh
bash tests/test_workflow.sh
BRAINOMICS_FULL_TESTS=1 bash tests/test_workflow.sh
python3 tests/test_package_contract.py
```

The default workflow check covers syntax, data/script consistency, interfaces,
and the available local package contract. `BRAINOMICS_FULL_TESTS=1` also runs
small synthetic tests of age mapping, LISI method identity, R/Python export and
readback, and PCA/RPCA/Harmony, neighbor and clustering algorithms. These tests
require the configured environment and exercise small datasets; they do not
rerun the full atlas or train scVI. Environment installation and verification
are separate from executing these tests. See
[Software checks](environment/REPRODUCIBILITY.md#software-checks) for details.

## Deposit package

The deposit package is stored outside the repository.
It holds 84 files (82 payload, about 22 GB): `expression/`, `embeddings/`, `validation/`, `metadata/`, `provenance/`, `scripts/`, `README.md`, `md5sum.txt`.

## Licence

Code is released under MIT (see `LICENSE`). GSE168408, Li_et_al_2018, Ma_et_al_2022, SomaMut, PRJCA015229 and ROSMAP were obtained as processed data and accompanying public annotations through author-provided download links or open repository records, without an additional project-specific access application. These public versions were harmonized and included with source attribution and documentation of their release-specific conditions. Actual routes, available license and policy information, and any uncertainty about the scope of terms are recorded in `data/source_access_summary.tsv` and `data/source_redistribution_evidence.tsv`. The code license does not relicense third-party data. Conditions for other files from the same study, such as raw sequencing or complete clinical metadata, are documented separately from those for the public processed inputs.
