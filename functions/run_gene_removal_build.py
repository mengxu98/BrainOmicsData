"""Sequential, resumable real matrix construction; uses existing installed R and sparse chunks."""
import csv
from decimal import Decimal
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess

run_value = os.environ.get('BRAINOMICS_RUN_ROOT')
assert run_value, 'BRAINOMICS_RUN_ROOT must name the run directory holding the source inputs and gene-panel policy'
run = Path(run_value).resolve()
out = run / 'run/gene_removal_build'
out.mkdir(parents=True, exist_ok=True)
lock = (out/'runner.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
resource.setrlimit(resource.RLIMIT_AS, (96 * 1024**3, 96 * 1024**3))
os.nice(19)
repo = Path(__file__).resolve().parents[1]
script = repo/'functions/build_clean_source_counts.R'
assert script.is_file(), f'Assembly source is missing: {script}'
list_file = Path(os.environ.get(
    'BRAINOMICS_DATASET_LIST', repo/'data/source_redistribution_evidence.tsv'))
assert list_file.is_file(), f'dataset list not found: {list_file}'
lines = list_file.read_text().splitlines()
datasets = lines if list_file.suffix == '.txt' else [line.split('	')[0] for line in lines[1:] if line]
assert len(datasets) == len(set(datasets)) == 22, datasets
assert len(datasets) == len(set(datasets)) == 22
env = dict(os.environ, BRAINOMICS_REPO_ROOT=str(repo), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
completed = []
state = dict(pid=os.getpid(), counts_modified=True, cell_membership_modified=False,
             sources=22, memory_limit_gib=96, sequential_sources=True)


def status(stage, **extra):
    state.update(state=stage, updated_utc=datetime.now(timezone.utc).isoformat(),
                 completed_sources=len(completed), **extra)
    temp = out/'status.json.tmp'
    temp.write_text(json.dumps(state, indent=2)+'\n'); temp.replace(out/'status.json')


try:
    for ds in datasets:
        destination = out/ds; destination.mkdir(exist_ok=True)
        marker = destination/'SUCCESS.json'
        status('BUILDING_CLEANED_COUNTS', current_source=ds)
        if not marker.exists():
            with (destination/'driver.log').open('a') as log:
                subprocess.run([os.environ.get('BRAINOMICS_RSCRIPT', 'Rscript'), str(script), str(run), ds],
                               cwd=repo, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        result = json.loads(marker.read_text())
        assert result['dataset'] == ds
        assert result['all_cells_visited_once'] and result['counts_modified']
        assert result['no_unsupported_positive_addition'] and not result['cell_membership_modified']
        completed.append(result)
    assert sum(r['reference_cells'] for r in completed) == 2602031
    status('SELECTING_FINAL_GENE_PANEL', current_source=None)
    maximum, total = {}, {}
    for ds in datasets:
        with (out/ds/'gene_detection.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        assert len(rows) == len({r['Gene_Key'] for r in rows}) == 34387
        for r in rows:
            value = Decimal(r['Detected_Cells'])
            assert value.is_finite() and value >= 0 and value == value.to_integral_value()
            key, n = r['Gene_Key'], int(value)
            maximum[key] = max(maximum.get(key, 0), n)
            total[key] = total.get(key, 0) + n
    with (run/'inputs/plan/selected_gene_panel_34387.tsv').open() as f:
        candidates = list(csv.DictReader(f, delimiter='\t'))
    fields = list(candidates[0]) + ['Cleaned_Max_Source_Detection', 'Cleaned_Total_Detection']
    retained, excluded = [], []
    for row in candidates:
        key = row['Gene_Key']
        row.update(Cleaned_Max_Source_Detection=maximum[key], Cleaned_Total_Detection=total[key])
        (retained if maximum[key] >= 10 else excluded).append(row)
    for name, rows in [('final_gene_panel.tsv', retained), ('genes_removed_after_count_cleaning.tsv', excluded)]:
        with (out/name).open('w') as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter='\t'); w.writeheader(); w.writerows(rows)
    summary = dict(reference_cells=2602031, sources=22, candidate_genes=34387,
                   final_genes=len(retained), removed_genes=len(excluded),
                   zero_count_cells_before_final_panel=sum(r['zero_count_cells'] for r in completed),
                   source_results=completed,
                   matrix_stage='Source count chunks complete; apply final panel to assembly and verify zero-count cells before integration.')
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    status('CLEANED_SOURCE_COUNTS_COMPLETE', current_source=None, final_genes=len(retained))
except Exception as error:
    status('FAILED', error=str(error))
    raise
