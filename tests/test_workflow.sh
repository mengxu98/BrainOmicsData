#!/usr/bin/env bash
# Consistency checks for the repository: syntax, data/script agreement and the
# released package contract.  Synthetic fixtures that exercise pipeline logic
# are only run with BRAINOMICS_FULL_TESTS=1.

set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"
source "$repo_dir/environment/activate.sh"

if [ -n "${BRAINOMICS_PYTHON:-}" ]; then
  python="$BRAINOMICS_PYTHON"
else
  python=python3
fi
export PYTHONDONTWRITEBYTECODE=1
test_dir="$(mktemp -d)"
trap 'rm -rf "$test_dir"' EXIT

# 1. Syntax: Python without py_compile, R parse, shell -n.
"$python" - "$test_dir/sources.txt" <<'PY'
import ast
from pathlib import Path
import shutil
import subprocess
import sys

repo = Path.cwd()
git_root = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True, text=True) if shutil.which('git') else None
if git_root is not None and git_root.returncode == 0 and Path(git_root.stdout.strip()).resolve() == repo.resolve():
    result = subprocess.run(['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'],
                            check=True, capture_output=True, text=True)
    paths = [Path(name) for name in result.stdout.split('\0') if name]
else:
    roots = ['processing', 'integration', 'tests', 'functions', 'environment', 'analysis',
             'annotation', 'plotting', 'sciencedb', 'download', 'hpc',
             'provenance/clustering_sensitivity/scripts']
    paths = list(repo.glob('*.sh')) + [path for name in roots for path in Path(name).rglob('*')]
sources = sorted({path for path in paths if path.is_file() and path.suffix in {'.py', '.R', '.sh', '.sbatch'}})
Path(sys.argv[1]).write_text(''.join(str(path) + '\n' for path in sources))
for path in sources:
    if path.suffix == '.py':
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
PY

"$BRAINOMICS_RSCRIPT" --vanilla - "$test_dir/sources.txt" <<'RS'
files <- readLines(commandArgs(TRUE)[1L])
files <- files[grepl("[.]R$", files)]
invisible(lapply(files, parse))
RS

while IFS= read -r script; do
  case "$script" in *.sh|*.sbatch) bash -n "$script" ;; esac
done < "$test_dir/sources.txt"

# 2. Data and script agreement: cohort, registry columns, download scripts,
#    package replay, driver references and environment locks.
"$python" tests/check_consistency.py

# 3. Download interface and the exact GSE168408 record.
for dataset in \
  AllenM1 \
  EGAS00001006537 \
  GSE168408 \
  Wang_2025 \
	GSE294786 \
	Velmeshev_2023; do
	[ -f "download/${dataset}.sh" ] || {
		echo "Missing download script for $dataset" >&2
		exit 1
	}
  download_records="$(bash "download/${dataset}.sh" --list)"
  if [ "$(printf '%s\n' "$download_records" | sed -n '1p')" != \
       "source_url|file_name|expected_bytes|sha256" ]; then
    echo "Invalid download-record interface for $dataset" >&2
    exit 1
  fi
done

source functions/utils.sh
gse168408_record="$(
  bash download/GSE168408.sh \
    --record RNA-all_full-counts-and-downsampled-CPM.h5ad
)"
expected_gse168408_record="https://storage.googleapis.com/neuro-dev/Processed_data/RNA-all_full-counts-and-downsampled-CPM.h5ad|RNA-all_full-counts-and-downsampled-CPM.h5ad|4963815084|5ca0b1ebe638f13e4302b721b80f41755a759c87f7fe43da6aeb7e8f9d0fd9ce"
if [ "$gse168408_record" != "$expected_gse168408_record" ]; then
  echo "GSE168408 processing input differs from its download record" >&2
  exit 1
fi

# 4. Script interfaces that failed silently before.
printf 'cross-platform-size\n' > "$test_dir/fixture.txt"
gzip -c "$test_dir/fixture.txt" > "$test_dir/fixture.txt.gz"
verify_file_integrity \
  "$test_dir/fixture.txt.gz" \
  "$(wc -c < "$test_dir/fixture.txt.gz" | tr -d ' ')"

for table in \
  data/donor_crosswalk.tsv \
  data/source_access_summary.tsv; do
  awk -F '\t' '
  NR == 1 { expected = NF }
  NF != expected {
    printf "%s:%d expected=%d got=%d\n", FILENAME, NR, expected, NF
    bad = 1
  }
  END { exit bad }
  ' "$table"
done

if HPC_PARALLEL_STREAMS=2 \
  bash hpc/pull_source_file.sh \
    EGAD00001006049 HumanFetalBrainPool.h5 >/dev/null 2>&1; then
  echo "Parallel transfer must require an explicit control-path prefix" >&2
  exit 1
fi

"$python" - "$repo_dir" "$test_dir" <<'PY'
import os
from pathlib import Path
import shlex
import subprocess
import sys
import importlib.util
from collections import defaultdict
import zipfile
import gzip
import json

repo, work = map(Path, sys.argv[1:])
stub = work/'sbatch'
stub.write_text('#!/bin/bash\nprintf "%s\\0" "$@" > "$BRAINOMICS_TEST_SUBMISSION"\n'
                'exit "${BRAINOMICS_TEST_SUBMIT_EXIT:-0}"\n')
stub.chmod(0o755)
record = work/'submission.bin'
env = dict(os.environ, PATH=str(work)+os.pathsep+os.environ['PATH'],
           BRAINOMICS_EXECUTOR='slurm', BRAINOMICS_HPC_PARTITION='cpu',
           BRAINOMICS_HPC_LOG_DIR=str(work/'logs'),
           BRAINOMICS_CPU_CPUS='4', BRAINOMICS_CPU_MEMORY='12G',
           BRAINOMICS_TEST_SUBMISSION=str(record))
cases = [('03_datasets_preprocessing.sh', ['T', 'GSE67835']),
         ('04_source_counts_assembly.sh', ['--check-inputs']),
         ('10_analysis_figures.sh', []), ('11_sciencedb_export.sh', [])]
for script, args in cases:
    result = subprocess.run(['bash', str(repo/script), *args], env=env,
                            capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, (script, result.stderr)
    submitted = record.read_bytes().decode().rstrip('\0').split('\0')
    assert '--wait' in submitted and '--export=ALL' in submitted
    assert '--cpus-per-task=4' in submitted and '--mem=12G' in submitted
    wrapped = next(a.split('=', 1)[1] for a in submitted if a.startswith('--wrap='))
    assert shlex.split(wrapped) == ['exec', 'env', 'BRAINOMICS_EXECUTOR=local',
                                    'bash', str(repo/script), *args]
failed = subprocess.run(['bash', str(repo/cases[0][0]), *cases[0][1]],
                        env=dict(env, BRAINOMICS_TEST_SUBMIT_EXIT='7'),
                        capture_output=True, text=True, timeout=20)
assert failed.returncode != 0
assert 'submission or execution failed' in failed.stderr
print('CPU stage submission, argument forwarding and failure propagation passed')
spec = importlib.util.spec_from_file_location('source_counts', repo/'functions/run_gene_removal_build.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
context = dict(identifiers=defaultdict(set, {'ENSG000001': {'HGNC:1'}, 'ENSG000002': {'HGNC:2'}}))
row = dict(Original_Feature_ID='old_name', Gene_Key='HGNC:1', Source_Ensembl='ENSG000001.5')
assert module.source_identity(row, 'GSE67835', context, {})[0]
assert not module.source_identity(dict(row, Source_Ensembl='ENSG000002'), 'GSE67835', context, {})[0]
assert not module.source_identity(dict(row, Source_Ensembl=''), 'GSE67835', context, {})[0]
assert module.source_identity(dict(row, Original_Feature_ID='ENSG000001.2', Source_Ensembl=''),
                              'GSE67835', context, {})[0]
assert not module.source_identity(dict(row, Source_Ensembl='ENSG000001;ENSG000002'),
                                  'GSE67835', context, {})[0]
reference = work/'Allen_GRCh38p2_reference.zip'
with zipfile.ZipFile(reference, 'w') as archive:
    archive.writestr('rsem_GRCh38.p2.gtf',
                     '1\tRefSeq\tgene\t1\t10\t.\t+\t.\tgene_id "123"; gene_symbol "old_name";\n')
context.update(reference_dir=work, reported={},
               entrez=defaultdict(set, {'123': {'HGNC:1'}}))
supported, ids = module.source_identity(row, 'AllenM1', context, {})
assert supported and ids == {'NCBI:123'}
print('Gene identity conflicts, missing IDs, embedded IDs and Allen source namespaces passed')
count_features = work/'source_counts.tsv'
count_features.write_text('Geneid\tcell\nENSG000001.5|supported_name\t3\n'
                          'ENSG000002|conflicting_name\t4\n'
                          'ENSG000001|merged_name\t1\nENSG000002|merged_name\t2\n')
source_ids = module.count_feature_ids(count_features)
source_row = dict(Original_Feature_ID='supported-name', Gene_Key='HGNC:1', Source_Ensembl='')
assert module.source_identity(source_row, 'GSE81475', context, source_ids)[0]
assert not module.source_identity(dict(source_row, Original_Feature_ID='conflicting-name'),
                                  'GSE81475', context, source_ids)[0]
assert not module.source_identity(dict(source_row, Original_Feature_ID='merged-name'),
                                  'GSE81475', context, source_ids)[0]
assert source_ids['supported-name'] == {'ENSG000001'}
print('Original count-row identities survive symbol aggregation; conflicting IDs remain excluded')
gene_refs = work/'gene_refs'
gene_refs.mkdir()
(gene_refs/'hgnc_complete_set.txt').write_text(
    'hgnc_id\tsymbol\tstatus\tprev_symbol\talias_symbol\tensembl_gene_id\n'
    'HGNC:1\tCURRENT\tApproved\t\tOLD_ALIAS\tENSG000001\n'
    'HGNC:2\tOTHER_CURRENT\tApproved\t\t\tENSG000005\n')
for filename, header in [('Homo_sapiens.gene_info.gz', '#tax_id\tGeneID\tdbXrefs\n'),
                          ('gene_history.gz', '#tax_id\tGeneID\tDiscontinued_GeneID\n')]:
    with gzip.open(gene_refs/filename, 'wt') as handle:
        handle.write(header)
for version in (17, 19, 24, 27, 32, 38, 43, 50):
    rows = {17: [('ENSG000002', 'OLD_ALIAS')], 19: [('ENSG000004', 'OTHER_CURRENT')],
            50: [('ENSG000001', 'CURRENT'), ('ENSG000003', 'CURRENT'),
                 ('ENSG000005', 'OTHER_CURRENT')]}.get(version, [])
    with gzip.open(gene_refs/f'gencode.v{version}.annotation.gtf.gz', 'wt') as handle:
        for identifier, name in rows:
            hgnc_tag = ' hgnc_id "HGNC:1";' if identifier == 'ENSG000001' else ''
            handle.write(f'1\tref\tgene\t1\t10\t.\t+\t.\tgene_id "{identifier}"; gene_name "{name}";{hgnc_tag}\n')
reference_context = module.gene_reference_context(gene_refs,
    [dict(Gene_Key='HGNC:1', Ensembl_IDs='ENSG000001')])
assert not reference_context['identifiers']['ENSG000002']
assert reference_context['identifiers']['ENSG000004'] == {'HGNC:2'}
assert reference_context['identifiers']['ENSG000003'] == {'ENSEMBL:ENSG000003'}
gencode_links = defaultdict(set)
module.gtf_gene_ids(gene_refs/'gencode.v50.annotation.gtf.gz', gencode_links)
assert gencode_links['ENSG000001'] == {'HGNC:1'}
assert not gencode_links['ENSG000003']
print('Historical source IDs require official reference names; aliases alone cannot resolve them')
source_run = work/'source_build'
rscript = os.environ['BRAINOMICS_RSCRIPT']
subprocess.run([rscript, '--vanilla', '-', str(source_run)], input='''
suppressPackageStartupMessages({library(SeuratObject); library(Matrix); library(data.table)})
run <- commandArgs(TRUE)[1L]
source <- file.path(run, 'inputs/processed/GSE67835')
plan <- file.path(run, 'inputs/plan')
dir.create(source, recursive=TRUE); dir.create(plan, recursive=TRUE)
counts <- Matrix(matrix(c(2,10,0, 0,10,3, 1,10,4), nrow=3), sparse=TRUE,
                 dimnames=list(c('CURRENT','CONFLICTING','OTHER-CURRENT'), c('a','b','c')))
saveRDS(CreateSeuratObject(counts), file.path(source, 'object.rds'))
fwrite(data.table(Format_Type='single_object', Full_Object='object.rds', Canonical_Metadata='metadata.tsv'),
       file.path(source, 'processed_object_format.tsv'), sep='\\t')
writeLines('Cell', file.path(source, 'metadata.tsv'))
genes <- data.table(Gene_Key=paste0('HGNC:', seq_len(34387)), Ensembl_IDs='')
genes$Ensembl_IDs[1:2] <- c('ENSG000001','ENSG000005')
fwrite(genes, file.path(plan, 'selected_gene_panel_34387.tsv'), sep='\\t')
fwrite(data.table(Dataset='GSE67835', Cells=paste0('GSE67835::', c('a','b','c'))),
       file.path(plan, 'retained_cells_exploration.tsv.gz'), sep='\\t')
dir.create(file.path(run, 'code'))
fwrite(data.table(Dataset='GSE67835', Original_Feature_ID=rownames(counts),
                 Existing_Feature=rownames(counts), Gene_Key=c('HGNC:1','HGNC:1','HGNC:2'),
                 Source_Ensembl=c('ENSG000001','ENSG000005','ENSG000005'), Action='candidate'),
       file.path(run, 'code/source_row_audit_map.tsv.gz'), sep='\\t')
fwrite(genes[1:24659], file.path(run, 'named_gene_panel_24659.tsv'), sep='\\t')
writeLines('{}', file.path(run, 'named_gene_panel_policy.json'))
''', text=True, check=True, capture_output=True)
source_env = dict(os.environ, BRAINOMICS_RUN_ROOT=str(source_run),
                  BRAINOMICS_EXECUTOR='local', BRAINOMICS_DATA_ROOT=str(work),
                  BRAINOMICS_GENE_REFERENCE_DIR=str(gene_refs), BRAINOMICS_REFERENCE_RUN_ROOT='')
checked = subprocess.run(['bash', str(repo/'04_source_counts_assembly.sh'), '--check-inputs',
                          '--dataset', 'GSE67835'], env=source_env, capture_output=True, text=True)
assert checked.returncode == 0, checked.stderr
rebuilt = subprocess.run([sys.executable, '-c',
                         "import resource, runpy, sys\n"
                         "if sys.platform == 'darwin': resource.setrlimit = lambda *args: None\n"
                         "sys.argv = sys.argv[1:]\n"
                         "runpy.run_path(sys.argv[0], run_name='__main__')",
                         str(repo/'functions/run_gene_removal_build.py'), '--dataset', 'GSE67835'],
                         env=source_env, capture_output=True, text=True)
assert rebuilt.returncode == 0, rebuilt.stderr + '\n'.join(
    path.read_text() for path in source_run.rglob('driver.log'))
policy = module.read_rows(source_run/'run/source_gene_policy/GSE67835/applied_gene_row_exclusions.tsv')
assert [row['Original_Feature_ID'] for row in policy] == ['CONFLICTING']
assert not (source_run/'inputs/input_manifest.json').exists()
subprocess.run([rscript, '--vanilla', '-', str(source_run)], input='''
suppressPackageStartupMessages(library(Matrix))
run <- commandArgs(TRUE)[1L]
counts <- readRDS(file.path(run, 'run/gene_removal_build/GSE67835/counts_00001.rds'))
stopifnot(identical(dim(counts), c(34387L,3L)),
          identical(colnames(counts), paste0('GSE67835::', c('a','b','c'))),
          all(as.matrix(counts[1:2,]) == matrix(c(2,0,0,3,1,4), nrow=2)),
          sum(counts[-c(1,2),]) == 0)
''', text=True, check=True, capture_output=True)
print('Source counts rebuild without hash manifests or a historical baseline; scientific exclusions and counts agree')
PY

mkdir -p "$test_dir/results/scvi_input"
ln -s /usr/bin/true "$test_dir/python3"
BRAINOMICS_SCVI_PYTHON="$test_dir/python3" BRAINOMICS_RESULTS_DIR="$test_dir/results" \
  bash 07_scvi.sh F >/dev/null

# 5. Released package contract (structure, hashes, table schemas, wording).
"$python" tests/test_package_contract.py

if find tests processing -type d -name __pycache__ -print -quit | grep -q .; then
  echo "Python cache directories must not be retained in the repository" >&2
  exit 1
fi

# 6. Optional synthetic fixtures (pipeline logic, not part of the default run).
if [ "${BRAINOMICS_FULL_TESTS:-0}" = "1" ]; then
  "$BRAINOMICS_RSCRIPT" --vanilla tests/test_celltype_metadata.R
  "$BRAINOMICS_RSCRIPT" --vanilla tests/test_integration_algorithms.R
  "$BRAINOMICS_RSCRIPT" --vanilla tests/test_lisi_neighbors.R
  "$BRAINOMICS_RSCRIPT" --vanilla tests/test_reference_knn_projection.R
  "$BRAINOMICS_RSCRIPT" --vanilla tests/test_query_normalization.R
  "$BRAINOMICS_RSCRIPT" --vanilla tests/test_age_interval_decisions.R
  "$BRAINOMICS_RSCRIPT" --vanilla tests/test_percent_mito_matrix.R
  if "$python" -c 'import h5py, scipy' >/dev/null 2>&1; then
    "$python" tests/test_percent_mito_h5ad.py
  else
    echo "percent_mito H5AD fixture skipped: h5py/scipy unavailable"
  fi
  if "$python" -c 'import h5py, numpy, pandas, scipy' >/dev/null 2>&1; then
    printf 'cell_id\tgene-a\tgene-b\ncell-1\t1\t0\ncell-2\t0\t2\ncell-3\t3\t4\n' \
      > "$test_dir/gse294786.tsv"
    "$python" processing/gse294786_tsv_to_h5ad.py \
      --counts "$test_dir/gse294786.tsv" \
      --output "$test_dir/gse294786.h5ad" \
      --audit-output "$test_dir/gse294786.audit.json" \
      --expected-cells 3 \
      --expected-features 2
    "$python" - "$test_dir/gse294786.h5ad" <<'PY'
import sys
import h5py

with h5py.File(sys.argv[1], "r") as handle:
    assert tuple(handle["X"].attrs["shape"]) == (3, 2)
    assert int(handle["X/indptr"][-1]) == 4
PY
  else
    echo "GSE294786 converter runtime skipped: h5py/pandas/scipy unavailable"
  fi
  "$BRAINOMICS_RSCRIPT" --vanilla tests/test_synthetic_pipeline.R
  echo "synthetic fixtures passed"
fi

echo "workflow checks passed"
