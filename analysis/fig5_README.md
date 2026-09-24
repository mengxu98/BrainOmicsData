# Figure 5 reproducible analysis

Run from the repository root. `analysis/fig5_rebuild.sh` extracts all 22
datasets from the prepared metadata, generates the four plotting tables,
then draws `figures/fig5a.pdf`, `fig5b.pdf`, `fig5c.pdf`, `fig5.pdf`, `fig5.svg`
and a 600 dpi `fig5.png`. Extraction needs the processed source objects and
the complete 22-source analysis inputs; it is not part of the regular
figure-only stage (`10_analysis_figures.sh`). The figure-only entry point is
`Rscript --vanilla plotting/fig5.R` once the four tables exist.

```sh
bash analysis/fig5_rebuild.sh
```

Optional paths:

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

`fig5_extract_full_panel.R` reads each source's raw RNA counts and canonical
feature crosswalk. It groups cells by source, canonical donor, age interval,
standardized brain region and adopted cell type, retaining unmeasured genes as
missing rather than zero. Every source cell must be covered exactly once.
Groups require at least 20 cells and 1,000 library counts. The 11 genes are
the prespecified reuse panel in the manuscript.

`fig5_prepare_sources.R` checks all 22 source-coverage files and 2,602,031 covered
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
also writes per-source coverage summaries and the preparer writes intermediate tables.
All source objects, generated tables and figures stay outside Git under the
repository's existing ignore policy; the analysis and plotting code are kept
in Git.
