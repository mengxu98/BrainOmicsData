# A human brain single-cell and single-nucleus transcriptomic resource across age intervals

BrainOmicsData is the code repository for this resource.

**2,602,031 cells or nuclei · 24,659 genes · 22 source datasets · 75 clusters · 12 cell types**

Version 2 is published at ScienceDB under the existing DOI [10.57760/sciencedb.41612](https://doi.org/10.57760/sciencedb.41612). The corresponding analysis code is identified by the [v2.0.0 tag](https://github.com/mengxu98/BrainOmicsData/tree/v2.0.0). Code maintenance does not replace the public ScienceDB files.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="figures/fig1-night.svg">
  <source media="(prefers-color-scheme: light)" srcset="figures/fig1.svg">
  <img alt="Overview of the integrated human brain atlas" src="figures/fig1.svg">
</picture>

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
| 11 | `11_analysis_figures.sh` | Figures from the retained analysis outputs |
| 12 | `12_sciencedb_upload.sh` | Transfer the sealed package |

`bash run_pipeline.sh` runs stages 01–12; `--list`, `--from`, `--to`, `--only` select subsets.
`BRAINOMICS_EXECUTOR=local` (default) runs stages in place; `BRAINOMICS_EXECUTOR=slurm` submits the matching `hpc/*.sbatch` job.

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

The manuscript figure entry points are `plotting/fig1.R` through
`plotting/fig5.R` and `plotting/figS1.R` through `plotting/figS4.R`. Each writes a
matching 600 dpi or higher PNG under `figures/`. Active panel PDFs use only
their panel identifiers (`fig1a.pdf`, `fig2f.pdf`, and so on); vector copies
and source tables remain in the same directory.
Figure 5's complete 22-source extraction, aggregation, source-level sign test,
fixed-context input, and panel assembly are documented in
[`analysis/fig5_README.md`](analysis/fig5_README.md). Run
`bash analysis/fig5_rebuild.sh` on the storage host to regenerate its source
tables and formal figure; stage 11 only redraws from existing source tables.
`plotting/` contains only figure-numbered
R scripts; shared configuration, panel assembly and other plotting support
scripts live under `functions/`. The S1 age-interval and S2 brain-region entry
points redraw their PDFs from the 2,602,031-cell metadata and RPCA embedding
before exporting PNGs.
The S3 cell-type entry point exports PNG from the retained PDF by default;
`Rscript plotting/figS3.R --redraw` draws it from the analysis inputs using
`functions/figS3_source.R`. The original drawing source is also available in
`analysis/revision_sources/plot_celltype_locations.R`.
S4 shows all 33 markers in the adopted list, with all 2,602,031 cells in each
FeatureDimPlot. Its complete drawing code is in `functions/figS4_source.R`;
`plotting/figS4.R` reads the full analysis inputs and exports the numbered PDF
and 600 dpi PNG. The plotting dependency is recorded in
`environment/scop-plotting.lock.tsv`, with source and license under
`environment/vendor/`. It uses the installed `scop` package or an explicitly
supplied `SCOP_SOURCE_PATH`; see
[`environment/REPRODUCIBILITY.md`](environment/REPRODUCIBILITY.md) for setup.

Local analysis inputs live under `results/`: `analysis_run/` (figure inputs), `frozen_run/` (stage-10 inputs), `annotation/`, `gene_reuse/`, `run_root/` (stage-04 run root; its large `inputs/`, `matrices/` and `run/` trees are kept on the storage host) and `work/` (stage-10 work directory). `results/`, `submission/`, `figures/` (except the light and dark Figure 1 SVGs) and the manuscript folders are outside version control. PNG exports remain local for manuscript preparation; GitHub displays the corresponding SVG for the reader's color theme.

The `integration_25` paths below are intermediate output locations for the numbered processing stages. The figures and annotation use the 2,602,031-cell inputs under `results/analysis_run/` and `results/annotation/`.

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
It holds 84 files (82 payload, about 22 GB): `expression/`, `embeddings/`, `validation/`, `metadata/`, `provenance/`, `scripts/`, `README.md`, `md5sum.txt`.

## Licence

MIT (see `LICENSE`).
