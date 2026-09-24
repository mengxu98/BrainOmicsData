# Clustering sensitivity analysis

This directory contains parameter records, overlap counts, per-cluster metrics,
and aggregate TSVs for five Louvain clustering settings on a fixed graph.
The analysis scripts are:

- `analysis/revision_sources/full_graph_stability.R`: executes the five settings.
- `analysis/revision_sources/summarize_stability.R`: summarizes partition and
  majority-matched cell-type stability.

## Settings and retained results

All five runs contain 2,602,031 cells and use the same graph MD5
(`f93f59d6fb1e79c0f413d7009491405a`) and reference-assignment MD5
(`5e1a06dc3744b5944125973295a6cc56`). Each run uses Louvain clustering
(`algorithm = 1`), `n.start = 10`, `n.iter = 10`, and
`group.singletons = TRUE`, without subsampling. The script runs `FindClusters()`
for the baseline before comparing its labels with the reference partition.

| Seed | Resolution | Clusters | ARI against the reference partition |
|---|---:|---:|---:|
| 20260730 | 2 | 75 | 1 |
| 42 | 2 | 77 | 0.754473754912277 |
| 2026 | 2 | 72 | 0.803025492214823 |
| 20260730 | 1.5 | 64 | 0.712367845263224 |
| 20260730 | 2.5 | 86 | 0.685946143287341 |

`stability/run_1/` through `stability/run_5/` contain the run parameters,
completion records, partition overlaps, and per-cluster summaries.
`stability/stability_summary.tsv`, `stability_by_dataset.tsv`,
`stability_by_celltype.tsv`, and `stability_by_cluster.tsv` contain the combined
results. `annotations/final_cluster_annotation.tsv` supplies the 75-cluster,
12-cell-type annotation used for majority matching.

Adjusted Rand index (ARI) compares partitions. Majority-matched cell-type
agreement measures stability conditional on the adopted annotation. These
analyses assess clustering seed and resolution on the same graph; integration
randomness, graph construction, and independent biological annotation accuracy
are outside their scope.

## Reproduction inputs and commands

The execution script accepts three arguments: an input base directory, an
existing output directory, and a task number (1–5). The input base must contain these files in the layout expected by the script:

```text
rpca_symmetric_weighted_knn_graph.rds
cluster_assignments.rds
```

Use a new output directory, preserve the original cell order, and use the
integration environment recorded in `environment/`. Run task 1 first; the
remaining tasks require its `BASELINE_REPLAY_PASS` marker.

```sh
mkdir -p /path/to/new-output
Rscript --vanilla analysis/revision_sources/full_graph_stability.R /path/to/input-base /path/to/new-output 1
# Repeat with task numbers 2, 3, 4 and 5 after task 1 passes.
cp provenance/clustering_sensitivity/annotations/final_cluster_annotation.tsv /path/to/new-output/
Rscript --vanilla analysis/revision_sources/summarize_stability.R /path/to/input-base /path/to/new-output
```

The full graph and per-cell assignments are required to rerun clustering and
summarization. For a smaller check, the following command recalculates ARI and
per-cluster metrics from the saved overlap counts and checks them against the
retained summaries. It does not rerun clustering or integration:

```sh
python3 - <<'PY'
from pathlib import Path
import runpy
root = Path('provenance/clustering_sensitivity')
verify = runpy.run_path(str(root / 'scripts/verify_stability_overlap.py'))['verify']
for task in range(1, 6):
    result = verify(root / 'stability' / f'run_{task}')
    print(task, result['state'], result['ARI'])
PY
```
