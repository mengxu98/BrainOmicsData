"""Assemble final source Seurat objects one at a time; no repeat filtering."""
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
out = run/'run/final_source_assembly'; out.mkdir(parents=True, exist_ok=True)
lock = (out/'runner.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
resource.setrlimit(resource.RLIMIT_AS, (96*1024**3, 96*1024**3)); os.nice(19)
repo = Path(__file__).resolve().parents[1]
script = repo/'functions/assemble_final_source.R'
assert script.is_file(), f'Assembly source is missing: {script}'
summary = json.loads((run/'run/gene_removal_build/summary.json').read_text())
assert summary['reference_cells'] == 2602031 and summary['sources'] == 22
list_file = Path(os.environ.get(
    'BRAINOMICS_DATASET_LIST', repo/'data/source_redistribution_evidence.tsv'))
assert list_file.is_file(), f'dataset list not found: {list_file}'
lines = list_file.read_text().splitlines()
datasets = lines if list_file.suffix == '.txt' else [line.split('	')[0] for line in lines[1:] if line]
assert len(datasets) == len(set(datasets)) == 22, datasets
env = dict(os.environ, BRAINOMICS_REPO_ROOT=str(repo), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
completed = []
state = dict(pid=os.getpid(), sources=22, genes=summary['final_genes'], memory_limit_gib=96)


def status(stage, **extra):
    state.update(state=stage, updated_utc=datetime.now(timezone.utc).isoformat(),
                 completed_sources=len(completed), **extra)
    temp=out/'status.json.tmp'; temp.write_text(json.dumps(state,indent=2)+'\n'); temp.replace(out/'status.json')


try:
    for ds in datasets:
        destination = run/'matrices/sources'/ds; destination.mkdir(parents=True,exist_ok=True)
        marker = destination/'SUCCESS.json'
        status('ASSEMBLING_FINAL_SOURCE_OBJECTS', current_source=ds)
        if not marker.exists():
            with (destination/'driver.log').open('a') as log:
                subprocess.run([os.environ.get('BRAINOMICS_RSCRIPT', 'Rscript'), str(script), str(run), ds],
                               cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        row=json.loads(marker.read_text())
        assert row['genes'] == summary['final_genes'] and row['zero_count_cells'] == 0
        assert not row['cell_membership_modified']
        completed.append(row)
    assert sum(r['reference_cells'] for r in completed) == 2602031
    (out/'summary.json').write_text(json.dumps(dict(reference_cells=2602031,genes=summary['final_genes'],
        sources=22,zero_count_cells=0,source_results=completed,
        next_stage='Combine accepted source objects in reference_datasets order, then matrix acceptance and integration.'),indent=2)+'\n')
    status('FINAL_SOURCE_OBJECTS_COMPLETE', current_source=None)
except Exception as error:
    status('FAILED', error=str(error)); raise
