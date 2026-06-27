#!/usr/bin/env bash
set -euo pipefail

source "functions/utils.sh"

check_command Rscript

overwrite="${1:-F}"

should_process() {
  local target_file="$1"
  if [[ "$overwrite" =~ ^([Tt]|[Tt][Rr][Uu][Ee]|1)$ ]]; then
    return 0
  fi
  if [ ! -f "$target_file" ]; then
    return 0
  fi
  return 1
}

if should_process "figures/datasets/group_heatmap_markergenes.pdf"; then
  log_message "Running integrated atlas plotting..."
  Rscript plotting/datasets.R
  log_message "Integrated atlas plotting completed successfully!" --message-type success
else
  log_message "Integrated atlas plotting outputs already exist!"
fi

if should_process "figures/lisi/lisi_plot_2.pdf"; then
  log_message "Running LISI plotting..."
  Rscript plotting/lisi_plot.R
  log_message "LISI plotting completed successfully!" --message-type success
else
  log_message "LISI plotting outputs already exist!"
fi
