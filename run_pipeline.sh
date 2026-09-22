#!/usr/bin/env bash
# Run the numbered deposit drivers in order.
#
# Usage:
#   bash run_pipeline.sh                    # stages 01 -> 12
#   bash run_pipeline.sh --from 05 --to 09  # a range
#   bash run_pipeline.sh --only 07,11       # selected stages
#   bash run_pipeline.sh --list             # stage numbers and files
#
# BRAINOMICS_EXECUTOR=local (default) runs every stage in place.
# BRAINOMICS_EXECUTOR=hpc submits hpc/*.sbatch and waits for each job.
#
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
cd "$BRAINOMICS_REPO_ROOT"

usage() {
  sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'
}

from=""
to=""
only=""
extra=()

while [ "$#" -gt 0 ]; do
  case "$1" in
    --from)
      [ "$#" -ge 2 ] || brainomics_die "--from needs a stage number"
      from="$2"
      shift 2
      ;;
    --to)
      [ "$#" -ge 2 ] || brainomics_die "--to needs a stage number"
      to="$2"
      shift 2
      ;;
    --only)
      [ "$#" -ge 2 ] || brainomics_die "--only needs stage numbers"
      only="$2"
      shift 2
      ;;
    --list)
      for script in [0-9][0-9]_*.sh; do
        printf '%s\t%s\n' "${script%%_*}" "$script"
      done
      exit 0
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      extra=("$@")
      break
      ;;
    *)
      brainomics_die "unknown argument: $1"
      ;;
  esac
done

stages=()
for script in [0-9][0-9]_*.sh; do
  number="${script%%_*}"
  if [ -n "$only" ]; then
    case ",$only," in
      *",$number,"*) ;;
      *) continue ;;
    esac
  else
    if [ -n "$from" ] && [ "$number" -lt "$from" ]; then
      continue
    fi
    if [ -n "$to" ] && [ "$number" -gt "$to" ]; then
      continue
    fi
  fi
  stages+=("$script")
done

[ "${#stages[@]}" -gt 0 ] || brainomics_die "no stage selected"

brainomics_log "executor=$BRAINOMICS_EXECUTOR stages=${#stages[@]}"
for script in "${stages[@]}"; do
  BRAINOMICS_STAGE="${script%.sh}"
  export BRAINOMICS_STAGE
  brainomics_log "starting $script"
  if [ "${#extra[@]}" -gt 0 ]; then
    bash "$script" "${extra[@]}"
  else
    bash "$script"
  fi
  brainomics_log "finished $script"
done
brainomics_log "pipeline completed"
