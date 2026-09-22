#!/usr/bin/env bash
# Consistency checks for the repository: syntax, data/script agreement and the
# released package contract.  Synthetic fixtures that exercise pipeline logic
# are only run with BRAINOMICS_FULL_TESTS=1.

set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"

if [ -n "${BRAINOMICS_PYTHON:-}" ]; then
  python="$BRAINOMICS_PYTHON"
elif [ -x "$repo_dir/.venv/bin/python" ]; then
  python="$repo_dir/.venv/bin/python"
else
  python=python3
fi
export PYTHONDONTWRITEBYTECODE=1

# 1. Syntax: Python without py_compile, R parse, shell -n.
"$python" - <<'PY'
import ast
from pathlib import Path

for root in (Path("processing"), Path("integration"), Path("tests")):
    for path in sorted(root.rglob("*.py")):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
PY

Rscript - <<'RS'
files <- sort(unique(c(
  list.files("functions", pattern = "[.]R$", full.names = TRUE),
  list.files("processing", pattern = "[.]R$", full.names = TRUE),
  list.files("integration", pattern = "[.]R$", full.names = TRUE),
  list.files("sciencedb", pattern = "[.]R$", full.names = TRUE),
  list.files("plotting", pattern = "[.]R$", full.names = TRUE),
  list.files("environment", pattern = "[.]R$", full.names = TRUE)
)))
invisible(lapply(files, parse))
RS

while IFS= read -r -d '' script; do
  bash -n "$script"
done < <(
  find . -type d \( -name .git -o -name .venv -o -name tmp \) -prune -o \
    -type f \( -name '*.sh' -o -name '*.sbatch' \) -print0
)

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
  [ -f "download/${dataset}.sh" ] || continue
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
test_dir="$(mktemp -d)"
trap 'rm -rf "$test_dir"' EXIT
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
  bash hpc/pull_source_file_from_hpc.sh \
    EGAD00001006049 HumanFetalBrainPool.h5 >/dev/null 2>&1; then
  echo "Parallel transfer must require an explicit control-path prefix" >&2
  exit 1
fi

mkdir -p "$test_dir/results/scvi_input"
ln -s /usr/bin/true "$test_dir/python3"
PATH="$test_dir:$PATH" BRAINOMICS_RESULTS_DIR="$test_dir/results" \
  bash 07_scvi.sh F >/dev/null

# 5. Released package contract (structure, hashes, table schemas, wording).
"$python" tests/test_package_contract.py

if find tests processing -type d -name __pycache__ -print -quit | grep -q .; then
  echo "Python cache directories must not be retained in the repository" >&2
  exit 1
fi

# 6. Optional synthetic fixtures (pipeline logic, not part of the default run).
if [ "${BRAINOMICS_FULL_TESTS:-0}" = "1" ]; then
  Rscript --vanilla tests/test_celltype_metadata.R
  Rscript --vanilla tests/test_annotation_inputs.R
  Rscript --vanilla tests/test_reference_knn_projection.R
  Rscript --vanilla tests/test_query_normalization.R
  Rscript --vanilla tests/test_age_interval_decisions.R
  Rscript --vanilla tests/test_percent_mito_matrix.R
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
  Rscript tests/test_clean_room_fixture.R
  echo "synthetic fixtures passed"
fi

echo "workflow checks passed"
