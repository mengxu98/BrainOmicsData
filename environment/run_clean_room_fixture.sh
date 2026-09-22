#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
profile="${1:-integration}"
case "$profile" in
  integration|preprocessing|local-validation) ;;
  *)
    echo "profile must be integration, preprocessing or local-validation" >&2
    exit 2
    ;;
esac

timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
if [ -n "${BRAINOMICS_REPRO_OUTPUT:-}" ]; then
  output_root="$BRAINOMICS_REPRO_OUTPUT"
elif [ "$profile" = "integration" ]; then
  output_root="$repo_dir/../../data/BrainOmicsData/integration_25/evaluation/clean_room"
else
  output_root="$repo_dir/results/clean_room"
fi
output_dir="$output_root/${profile}_${timestamp}"
if [ -e "$output_dir" ]; then
  echo "clean-room output already exists: $output_dir" >&2
  exit 2
fi
mkdir -p "$output_dir"
log_file="$output_dir/execution.log"
success_file="$output_dir/_SUCCESS"

exec > >(tee "$log_file") 2>&1
echo "Status=running"
echo "Profile=$profile"
echo "Started_UTC=$timestamp"
echo "Host=$(hostname)"
echo "Repository=$repo_dir"
if (cd "$repo_dir" && git rev-parse --is-inside-work-tree >/dev/null 2>&1); then
  git_head="$(cd "$repo_dir" && git rev-parse HEAD)"
  git_status="$(cd "$repo_dir" && git status --porcelain)"
  echo "Git_HEAD=$git_head"
  if [ -n "$git_status" ]; then
    echo "Git_Worktree=modified"
  else
    echo "Git_Worktree=clean"
  fi
else
  echo "Git_HEAD=not_available"
  echo "Git_Worktree=not_a_git_checkout"
fi

cd "$repo_dir"
export BRAINOMICS_REQUIRE_PYTHON_READER=true
if [ "${BRAINOMICS_REQUIRE_CLEAN_CHECKOUT:-false}" = "true" ]; then
  test "${git_status-unavailable}" = "" || { echo "Clean code snapshot required"; exit 1; }
fi
Rscript --vanilla environment/verify_environment.R "$repo_dir" "$profile"
if [ "${BRAINOMICS_VERIFY_SCVI_PYTHON:-false}" = "true" ]; then
  python_bin="python3"
  python_lock="environment/python-scvi-freeze.txt"
  "$python_bin" environment/verify_python_environment.py \
    --lock "$python_lock"
fi

bash tests/test_workflow.sh

find \
  functions processing integration annotation plotting sciencedb environment \
  tests hpc \
  -type f \
  \( -name '*.R' -o -name '*.py' -o -name '*.sh' -o -name '*.sbatch' \
     -o -name '*.tsv' -o -name '*.txt' -o -name '*.md' \) \
  -print | LC_ALL=C sort | while IFS= read -r file; do
    if command -v sha256sum >/dev/null 2>&1; then
      sha256sum "$file"
    else
      shasum -a 256 "$file"
    fi
  done > "$output_dir/code_sha256.tsv"

finished="$(date -u +%Y%m%dT%H%M%SZ)"
printf '%s\n' \
  "Status=complete" \
  "Profile=$profile" \
  "Started_UTC=$timestamp" \
  "Finished_UTC=$finished" \
  "Execution_Log=$(basename "$log_file")" \
  "Code_Manifest=code_sha256.tsv" > "$success_file"
echo "Status=complete"
echo "Finished_UTC=$finished"
