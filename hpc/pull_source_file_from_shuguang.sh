#!/usr/bin/env bash

set -euo pipefail

dataset="${1:?Usage: pull_source_file_from_hpc.sh <dataset> <filename> [expected_bytes]}"
filename="${2:?Usage: pull_source_file_from_hpc.sh <dataset> <filename> [expected_bytes]}"
requested_bytes="${3:-}"

case "$dataset" in
  AllenM1|EGAD00001006049|EGAS00001006537|GSE144136|GSE178175|\
  GSE202210|Catching_2026|Wang_2025|GSE294786|Velmeshev_2023)
    ;;
  *)
    echo "Unsupported dataset: $dataset" >&2
    exit 2
    ;;
esac
case "$filename" in
  ""|*/*|.|..|*[!A-Za-z0-9._-]*)
    echo "filename must be one safe basename without path traversal" >&2
    exit 2
    ;;
esac

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$repo_dir/functions/utils.sh"
data_root="${BRAINOMICS_DATA_ROOT:-$HOME/data/BrainOmicsData}"
hpc_data_root="/path/to/hpc/home/data/BrainOmicsData"
hpc_host="${HPC_HOST:-user@login-host}"
hpc_port="${HPC_PORT:-22}"
known_hosts="${HPC_KNOWN_HOSTS:-$HOME/.ssh/known_hosts_brainomics_hpc}"
control_path="${HPC_CONTROL_PATH:-}"
control_path_prefix="${HPC_CONTROL_PATH_PREFIX:-}"
parallel_streams="${HPC_PARALLEL_STREAMS:-1}"
audit_dir="$data_root/transfer_audit/additional_datasets"

case "$parallel_streams" in
  ""|*[!0-9]*)
    echo "HPC_PARALLEL_STREAMS must be an integer from 1 to 16" >&2
    exit 2
    ;;
esac
if [ "$parallel_streams" -lt 1 ] || [ "$parallel_streams" -gt 16 ]; then
  echo "HPC_PARALLEL_STREAMS must be an integer from 1 to 16" >&2
  exit 2
fi
if [ "$parallel_streams" -gt 1 ] && [ -z "$control_path_prefix" ]; then
  echo "HPC_CONTROL_PATH_PREFIX is required for parallel transfer" >&2
  exit 2
fi

download_record="$(
  dataset_download_record "$repo_dir" "$dataset" "$filename"
)" || exit 1
IFS='|' read -r \
  source_url record_filename expected_bytes expected_sha256 \
  <<< "$download_record"
if [ -n "$requested_bytes" ] && [ "$requested_bytes" != "$expected_bytes" ]; then
  echo \
    "Requested byte size differs from the download record for $dataset/$filename" \
    >&2
  exit 2
fi

ssh_options=(
  -o BatchMode=yes
  -o ServerAliveInterval=10
  -o ServerAliveCountMax=120
  -o StrictHostKeyChecking=yes
  -o "UserKnownHostsFile=$known_hosts"
  -p "$hpc_port"
)
if [ -n "$control_path" ]; then
  ssh_options+=(
    -o "ControlPath=$control_path"
    -o ControlMaster=no
  )
fi
source_file="$hpc_data_root/raw/$dataset/$filename"
destination_dir="$data_root/raw/$dataset"
destination_file="$destination_dir/$filename"
transfer_lock="${destination_file}.transfer.lock"

mkdir -p "$destination_dir" "$audit_dir"
if ! mkdir "$transfer_lock" 2>/dev/null; then
  echo "Another transfer already owns $dataset/$filename: $transfer_lock" >&2
  exit 75
fi
cleanup_transfer_lock() {
  rmdir "$transfer_lock" 2>/dev/null || true
}
trap cleanup_transfer_lock EXIT
trap 'exit 130' HUP INT TERM

source_bytes="$(
  ssh "${ssh_options[@]}" "$hpc_host" \
    "stat -c %s -- '$source_file'"
)"
if [ "$source_bytes" != "$expected_bytes" ]; then
  echo "Source byte-size mismatch for $dataset/$filename: expected $expected_bytes, observed $source_bytes" >&2
  exit 1
fi

parallel_range_transfer() {
  local prefix_bytes prefix_source_sha prefix_destination_sha
  local remaining_bytes chunk_bytes block_count index part
  local part_dir task_dir stage_file backup_file pid status
  local -a pids=()

  prefix_bytes=0
  if [ -f "$destination_file" ]; then
    prefix_bytes="$(stat -c %s "$destination_file")"
  fi
  if [ "$prefix_bytes" -gt "$expected_bytes" ]; then
    echo "Existing destination is larger than the source" >&2
    return 1
  fi
  if [ "$prefix_bytes" -gt 0 ]; then
    prefix_source_sha="$(
      ssh "${ssh_options[@]}" "$hpc_host" \
        "head -c '$prefix_bytes' -- '$source_file' | sha256sum" |
        awk '{print $1}'
    )"
    prefix_destination_sha="$(
      head -c "$prefix_bytes" "$destination_file" | sha256sum |
        awk '{print $1}'
    )"
    if [ "$prefix_source_sha" != "$prefix_destination_sha" ]; then
      backup_file="${destination_file}.invalid-prefix.$(date -u +%Y%m%dT%H%M%SZ)"
      mv "$destination_file" "$backup_file"
      prefix_bytes=0
    fi
  fi

  remaining_bytes=$((expected_bytes - prefix_bytes))
  if [ "$remaining_bytes" -eq 0 ]; then
    return 0
  fi
  chunk_bytes="${HPC_PARALLEL_CHUNK_BYTES:-67108864}"
  case "$chunk_bytes" in
    ""|*[!0-9]*)
      echo "HPC_PARALLEL_CHUNK_BYTES must be a positive integer" >&2
      return 2
      ;;
  esac
  if [ "$chunk_bytes" -lt 1048576 ]; then
    echo "HPC_PARALLEL_CHUNK_BYTES must be at least 1048576" >&2
    return 2
  fi
  block_count=$(((remaining_bytes + chunk_bytes - 1) / chunk_bytes))
  part_dir="${destination_file}.parallel-blocks"
  task_dir="$part_dir/tasks"
  mkdir -p "$part_dir" "$task_dir"

  # Each task is a relatively small, ordered byte range. Workers claim the
  # next available task atomically, so one slow WAN connection can delay only
  # one block instead of owning one eighth of a multi-gigabyte file.
  for ((index = 0; index < block_count; index++)); do
    printf '%d\n' "$index" > "$task_dir/todo.$(printf '%06d' "$index")"
  done

  parallel_range_worker() {
    local worker_index="$1"
    local socket task claim block block_name block_offset block_length
    local block_file block_bytes range_offset range_length source_part_sha
    local destination_part_sha log_file

    socket="${control_path_prefix}.${worker_index}"
    while :; do
      task="$(
        find "$task_dir" -maxdepth 1 -type f \
          -name 'todo.[0-9][0-9][0-9][0-9][0-9][0-9]' -print -quit
      )"
      [ -n "$task" ] || break
      claim="${task}.running.${worker_index}"
      mv "$task" "$claim" 2>/dev/null || continue
      block="$(cat "$claim")"
      block_name="$(printf '%06d' "$block")"
      block_offset=$((prefix_bytes + block * chunk_bytes))
      block_length=$chunk_bytes
      if [ $((block_offset + block_length)) -gt "$expected_bytes" ]; then
        block_length=$((expected_bytes - block_offset))
      fi
      block_file="$part_dir/part.$block_name"
      log_file="$part_dir/part.$block_name.log"
      block_bytes=0
      if [ -f "$block_file" ]; then
        block_bytes="$(stat -c %s "$block_file")"
      fi
      if [ "$block_bytes" -gt "$block_length" ]; then
        mv "$block_file" \
          "${block_file}.oversize.$(date -u +%Y%m%dT%H%M%SZ)"
        block_bytes=0
      fi
      if [ "$block_bytes" -gt 0 ]; then
        source_part_sha="$(
          ssh \
            -o BatchMode=yes \
            -o ServerAliveInterval=10 \
            -o ServerAliveCountMax=120 \
            -o StrictHostKeyChecking=yes \
            -o "UserKnownHostsFile=$known_hosts" \
            -o "ControlPath=$socket" \
            -o ControlMaster=no \
            -p "$hpc_port" \
            "$hpc_host" \
            "dd if='$source_file' bs=8388608 iflag=skip_bytes,count_bytes skip='$block_offset' count='$block_bytes' status=none | sha256sum" |
            awk '{print $1}'
        )"
        destination_part_sha="$(sha256sum "$block_file" | awk '{print $1}')"
        if [ "$source_part_sha" != "$destination_part_sha" ]; then
          mv "$block_file" \
            "${block_file}.invalid.$(date -u +%Y%m%dT%H%M%SZ)"
          block_bytes=0
        fi
      fi
      range_offset=$((block_offset + block_bytes))
      range_length=$((block_length - block_bytes))
      if [ "$range_length" -gt 0 ]; then
        ssh \
          -o BatchMode=yes \
          -o ServerAliveInterval=10 \
          -o ServerAliveCountMax=120 \
          -o StrictHostKeyChecking=yes \
          -o "UserKnownHostsFile=$known_hosts" \
          -o "ControlPath=$socket" \
          -o ControlMaster=no \
          -p "$hpc_port" \
          "$hpc_host" \
          "dd if='$source_file' bs=8388608 iflag=skip_bytes,count_bytes skip='$range_offset' count='$range_length' status=none" \
          >> "$block_file" 2> "$log_file"
      fi
      if [ "$(stat -c %s "$block_file")" != "$block_length" ]; then
        echo "Parallel block size mismatch: $block_file" >&2
        return 1
      fi
      mv "$claim" "$task_dir/done.$block_name"
    done
  }

  for ((index = 0; index < parallel_streams; index++)); do
    parallel_range_worker "$index" &
    pids+=("$!")
  done

  status=0
  for pid in "${pids[@]}"; do
    wait "$pid" || status=1
  done
  if [ "$status" -ne 0 ]; then
    echo "One or more parallel transfer ranges failed" >&2
    return 1
  fi

  stage_file="${destination_file}.parallel-stage"
  if [ "$prefix_bytes" -gt 0 ]; then
    cp --reflink=auto "$destination_file" "$stage_file"
  else
    truncate -s 0 "$stage_file"
  fi
  for ((index = 0; index < block_count; index++)); do
    part="$part_dir/part.$(printf '%06d' "$index")"
    if [ ! -f "$part" ]; then
      echo "Missing parallel block: $part" >&2
      return 1
    fi
    dd if="$part" of="$stage_file" bs=8388608 oflag=append conv=notrunc status=none
  done
  if [ "$(stat -c %s "$stage_file")" != "$expected_bytes" ]; then
    echo "Parallel stage byte-size mismatch" >&2
    return 1
  fi
  if [ "$(sha256sum "$stage_file" | awk '{print $1}')" != "$expected_sha256" ]; then
    echo "Parallel stage SHA-256 differs from the download record" >&2
    return 1
  fi
  if [ -f "$destination_file" ]; then
    backup_file="${destination_file}.partial.$(date -u +%Y%m%dT%H%M%SZ)"
    mv "$destination_file" "$backup_file"
  fi
  mv "$stage_file" "$destination_file"
  find "$part_dir" -maxdepth 1 -type f \
    \( -name 'part.[0-9][0-9][0-9][0-9][0-9][0-9]' -o \
       -name 'part.[0-9][0-9][0-9][0-9][0-9][0-9].log' \) -delete
  find "$task_dir" -maxdepth 1 -type f -name 'done.*' -delete
  rmdir "$task_dir" 2>/dev/null || true
  rmdir "$part_dir" 2>/dev/null || true
}

rsync_options=(
  -rlt
  --partial
  --append-verify
  --itemize-changes
  --stats
)
# Large plain-text matrices compress efficiently in transit. HDF5, gzip, zip,
# H5AD and RDS inputs remain byte-for-byte transfers because live throughput
# tests showed that extra HDF5 compression is CPU-bound on the source host.
case "$filename" in
  *.csv|*.tsv|*.txt|*.json)
    rsync_options+=(--compress)
    ;;
esac

if [ "$parallel_streams" -gt 1 ] && \
   { [ ! -f "$destination_file" ] || \
     [ "$(stat -c %s "$destination_file")" != "$expected_bytes" ]; }; then
  parallel_range_transfer
else
  rsync \
    "${rsync_options[@]}" \
    -e "ssh ${ssh_options[*]}" \
    "$hpc_host:$source_file" \
    "$destination_file"
fi

destination_bytes="$(stat -c %s "$destination_file")"
if [ "$destination_bytes" != "$expected_bytes" ]; then
  echo "Destination byte-size mismatch for $dataset/$filename: expected $expected_bytes, observed $destination_bytes" >&2
  exit 1
fi

source_sha256="$(
  ssh "${ssh_options[@]}" "$hpc_host" \
    "sha256sum -- '$source_file'" |
    awk '{print $1}'
)"
destination_sha256="$(sha256sum "$destination_file" | awk '{print $1}')"
if [ "$source_sha256" != "$destination_sha256" ]; then
  echo "SHA-256 mismatch for $dataset/$filename" >&2
  exit 1
fi
if [ "$source_sha256" != "$expected_sha256" ]; then
  echo "SHA-256 differs from the download record for $dataset/$filename" >&2
  exit 1
fi

# A failed direct download on 131 can leave these exact resumable companions.
# They are obsolete only after the complete destination has passed both the
# download-record and cross-host SHA-256 gates above.
rm -f "${destination_file}.tmp" "${destination_file}.resume"

safe_filename="${filename//[^A-Za-z0-9._-]/_}"
audit_file="$audit_dir/${dataset}__${safe_filename}.tsv"
temporary_audit="${audit_file}.tmp.$$"
{
  printf 'Dataset\tFilename\tExpected_Bytes\tSource_Bytes\tDestination_Bytes\tSource_SHA256\tDestination_SHA256\tVerified_UTC\tSource_Deleted\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tfalse\n' \
    "$dataset" \
    "$filename" \
    "$expected_bytes" \
    "$source_bytes" \
    "$destination_bytes" \
    "$source_sha256" \
    "$destination_sha256" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$temporary_audit"
mv "$temporary_audit" "$audit_file"

printf 'Verified complete transfer: %s/%s (%s bytes; SHA-256 %s)\n' \
  "$dataset" "$filename" "$expected_bytes" "$source_sha256"
