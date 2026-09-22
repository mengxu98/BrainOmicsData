# BrainOmicsData

Human brain single-cell and single-nucleus transcriptome resource.

**2,602,031 cells or nuclei · 24,659 genes · 22 source datasets · 75 clusters · 12 cell types**

Version 1 is published at ScienceDB (DOI [10.57760/sciencedb.41612](https://doi.org/10.57760/sciencedb.41612)); this repository holds the code for version 2.

![Overview of the integrated human brain atlas](figures/fig1.svg)

## Pipeline

Twelve numbered drivers run in order from the repository root:

| Stage | Driver | What it does |
|---|---|---|
| 01 | `01_environment.sh` | Verify the locked R and Python environments |
| 02 | `02_datasets_download.sh` | Download each source dataset |
| 03 | `03_datasets_preprocessing.sh` | Reconstruct and standardize each source |
| 04 | `04_source_counts_assembly.sh` | Build cleaned counts, apply the gene panel, accept the common matrices |
| 05 | `05_hvg_scvi_input.sh` | Select 3,000 HVGs and export the scVI input |
| 06 | `06_pca_harmony.sh` | Raw PCA, raw UMAP, Harmony |
| 07 | `07_scvi.sh` | scVI latent space |
| 08 | `08_rpca.sh` | RPCA latent space |
| 09 | `09_final_assembly.sh` | Assemble the integrated object and evaluate latent spaces |
| 10 | `10_sciencedb_export.sh` | Build the cell-level QC and annotation inputs when absent, then export and verify the deposit package |
| 11 | `11_analysis_figures.sh` | Annotation summaries, manuscript figures and tables |
| 12 | `12_sciencedb_upload.sh` | Transfer the sealed package |

`bash run_pipeline.sh` runs stages 01–12; `--list`, `--from`, `--to`, `--only` select subsets.
`BRAINOMICS_EXECUTOR=local` (default) runs stages in place; `BRAINOMICS_EXECUTOR=hpc` submits the matching `hpc/*.sbatch` job.

## Layout

| Path | Contents |
|---|---|
| `functions/` | Shared helpers, assembly runners, package exporters, normalization scripts |
| `integration/` | Integration and evaluation modules |
| `processing/`, `download/` | Per-dataset reconstruction and download scripts |
| `annotation/` | Annotation inputs and adopted-label export |
| `analysis/`, `plotting/` | Analysis summaries and figures |
| `sciencedb/` | Package metadata, manifest, reader and anonymization utilities |
| `results/` | Local analysis inputs used by the figure scripts (outside version control) |
| `environment/` | Package locks and restore/verification scripts |
| `hpc/` | HPC submission scripts (`hpc/*.sbatch`) and transfer helpers |
| `data/` | Source access table, donor crosswalks, feature metadata |
| `tests/` | Workflow tests and package contract checks |

Local analysis inputs live under `results/`: `analysis_run/` (figure inputs), `frozen_run/` (stage-10 inputs), `annotation/`, `gene_reuse/`, `run_root/` (stage-04 run root; its large `inputs/`, `matrices/` and `run/` trees are kept on the storage host) and `work/` (stage-10 work directory). `results/`, `submission/`, `figures/` (except the Figure 1 overview files) and the manuscript folders are outside version control.

## Environment variables

| Variable | Purpose | Used by |
|---|---|---|
| `BRAINOMICS_DATA_ROOT` | Dataset root; defaults to `<repository>/../../data/BrainOmicsData` | 03, 04, 10 |
| `BRAINOMICS_RESULTS_DIR` | Integration results; defaults to `$BRAINOMICS_DATA_ROOT/integration_25` | 05, 07, 09, 10 |
| `BRAINOMICS_RUN_ROOT` | Run directory holding `pipeline/` and `integration/` | 04-09 |
| `BRAINOMICS_FROZEN_DIR` | Frozen results of the export run | 10 |
| `BRAINOMICS_ANALYSIS_DIR` | Analysis outputs of the export run | 10 |
| `BRAINOMICS_PACKAGE_DIR` | Deposit package directory; defaults to `$BRAINOMICS_DATA_ROOT/ScienceDB` | 10 |
| `BRAINOMICS_WORK_DIR` | Export work directory (cell-level QC and annotation tables) | 10 |
| `BRAINOMICS_LISI_DIR` | Per-cell iLISI checkpoints | 10 |

Stage 10 also reads the frozen results and analysis outputs of the run that produced the release.

## Tests

```sh
bash tests/test_workflow.sh                           # syntax, data/script consistency, package contract
BRAINOMICS_FULL_TESTS=1 bash tests/test_workflow.sh   # add the synthetic fixtures (unit tests and clean-room end-to-end)
python3 tests/test_package_contract.py                # package contract only (add --strict for a full copy)
```

`tests/check_consistency.py` is the core of the default run: it compares the 22-dataset cohort
across the evidence table, `functions/integration.R`, the preprocessing driver, the HPC
dispatcher, the ScienceDB source-input contract and `processing/<dataset>.R`; checks the registry
columns and download scripts; replays `provenance/dataset_manifest.tsv` from the repository tables
and compares it with the released package; and verifies that every path referenced by the numbered
drivers exists.

## Deposit package

The sealed deposit package is kept on the storage host outside the repository.
It holds 85 files (83 payload, about 22 GB): `expression/`, `embeddings/`, `validation/`, `metadata/`, `provenance/`, `scripts/`, `README.md`, `md5sum.txt`.

## Licence

MIT (see `LICENSE`).
