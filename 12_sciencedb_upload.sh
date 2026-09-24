#!/usr/bin/env bash
# Transfer the sealed deposit package to the deposit server over FTP.
# Requires credentials on the command line; nothing is uploaded by default.
#
# Requires the sealed package, deposit-server credentials and network access.

set -euo pipefail

BRAINOMICS_STAGE=12_sciencedb_upload

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/functions/pipeline_lib.sh"
source "$SCRIPT_DIR/functions/log_message.sh"

usage() {
  cat <<'EOF'
Usage:
  bash 12_sciencedb_upload.sh --host HOST --port PORT --user USER --password PASSWORD [options]

Options:
  --local-dir DIR      Local package directory to upload.
                       Default: $BRAINOMICS_RUN_ROOT/package (or results/run_root/package)
  --remote-dir DIR     Remote FTP directory.
                       Default: human_brain_age_interval_sc_snRNAseq_dataset
  --host HOST          FTP host, for example ftp-upload.scidb.cn
  --port PORT          FTP port, for example 2121
  --user USER          FTP user name
  --password PASSWORD  FTP password
  --force              Upload even when the remote file has the same size
  --changed-by-md5     Download remote md5sum.txt and upload only files whose
                       local md5 differs from the remote manifest. md5sum.txt
                       is uploaded last when changed.
  --delete-remote-extra
                       Delete known obsolete files from previous package layouts
                       and remove the empty remote figures directory if present.
  -h, --help           Show this help

The FTP host, port, user and password can also be supplied as environment
variables: FTP_HOST, FTP_PORT, FTP_USER and FTP_PASSWORD. Credentials are not
written to files by this script.
EOF
}

LOCAL_DIR="${BRAINOMICS_PACKAGE_DIR:-${BRAINOMICS_RUN_ROOT:-$SCRIPT_DIR/results/run_root}/package}"
REMOTE_DIR="human_brain_age_interval_sc_snRNAseq_dataset"
FTP_HOST="${FTP_HOST:-}"
FTP_PORT="${FTP_PORT:-}"
FTP_USER="${FTP_USER:-}"
FTP_PASSWORD="${FTP_PASSWORD:-}"
FORCE=0
CHANGED_BY_MD5=0
DELETE_REMOTE_EXTRA=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --local-dir)
      LOCAL_DIR="$2"
      shift 2
      ;;
    --remote-dir)
      REMOTE_DIR="$2"
      shift 2
      ;;
    --host)
      FTP_HOST="$2"
      shift 2
      ;;
    --port)
      FTP_PORT="$2"
      shift 2
      ;;
    --user)
      FTP_USER="$2"
      shift 2
      ;;
    --password)
      FTP_PASSWORD="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --changed-by-md5)
      CHANGED_BY_MD5=1
      shift
      ;;
    --delete-remote-extra)
      DELETE_REMOTE_EXTRA=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      log_message "Unknown argument: {.arg $1}" --message-type error
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$FTP_HOST" || -z "$FTP_PORT" || -z "$FTP_USER" || -z "$FTP_PASSWORD" ]]; then
  log_message "Missing FTP credentials or endpoint." --message-type error
  usage >&2
  exit 2
fi

if [[ ! -d "$LOCAL_DIR" ]]; then
  log_message "Local directory does not exist: {.path $LOCAL_DIR}" --message-type error
  exit 1
fi
for required in md5sum.txt provenance/file_manifest.tsv; do
  if [[ ! -f "$LOCAL_DIR/$required" ]]; then
    log_message "Local package is missing {.file $required}." --message-type error
    exit 1
  fi
done

if ! command -v curl >/dev/null 2>&1; then
  log_message "{.pkg curl} is required." --message-type error
  exit 1
fi

LOCAL_DIR="$(cd "$LOCAL_DIR" && pwd)"
REMOTE_DIR="${REMOTE_DIR#/}"
REMOTE_DIR="${REMOTE_DIR%/}"

remote_size() {
  local rel="$1"
  local url="ftp://${FTP_HOST}:${FTP_PORT}/${REMOTE_DIR}/${rel}"
  curl --silent --show-error --fail --head \
    --connect-timeout 60 \
    --max-time 90 \
    --user "${FTP_USER}:${FTP_PASSWORD}" \
    "$url" 2>/dev/null |
    awk 'BEGIN{IGNORECASE=1} /^Content-Length:/ {gsub("\r", "", $2); print $2; exit}'
}

remote_files() {
  curl --silent --show-error --fail --list-only \
    --connect-timeout 60 \
    --user "${FTP_USER}:${FTP_PASSWORD}" \
    "ftp://${FTP_HOST}:${FTP_PORT}/${REMOTE_DIR}/" 2>/dev/null || true
}

delete_remote_file() {
  local rel="$1"
  curl --silent --show-error \
    --connect-timeout 60 \
    --user "${FTP_USER}:${FTP_PASSWORD}" \
    --quote "DELE /${REMOTE_DIR}/${rel}" \
    "ftp://${FTP_HOST}:${FTP_PORT}/" >/dev/null || true
}

delete_remote_dir() {
  local rel="$1"
  curl --silent --show-error \
    --connect-timeout 60 \
    --user "${FTP_USER}:${FTP_PASSWORD}" \
    --quote "RMD /${REMOTE_DIR}/${rel}" \
    "ftp://${FTP_HOST}:${FTP_PORT}/" >/dev/null || true
}

upload_file() {
  local rel="$1"
  local file="${LOCAL_DIR}/${rel}"
  local local_size remote_size_value url

  local_size="$(wc -c < "$file" | tr -d ' ')"
  remote_size_value="$(remote_size "$rel" || true)"

  if [[ "$FORCE" -eq 0 && -n "$remote_size_value" && "$remote_size_value" == "$local_size" ]]; then
    log_message "Skipping {.file $rel} ({.val {$local_size}} bytes)"
    return 0
  fi

  log_message "Uploading {.file $rel} local={.val {$local_size}} remote={.val {${remote_size_value:-NA}}}" --message-type running
  url="ftp://${FTP_HOST}:${FTP_PORT}/${REMOTE_DIR}/${rel}"
  curl \
    --fail \
    --ftp-create-dirs \
    --retry 100 \
    --retry-delay 60 \
    --connect-timeout 120 \
    --speed-time 900 \
    --speed-limit 256 \
    --user "${FTP_USER}:${FTP_PASSWORD}" \
    --upload-file "$file" \
    "$url"

  # The deposit server does not always answer SIZE immediately after a large
  # transfer. Retry a few times; treat a persistently unavailable size as a
  # warning instead of failing the whole run.
  local attempt=1
  while [[ "$attempt" -le 4 ]]; do
    remote_size_value="$(remote_size "$rel" || true)"
    if [[ -n "$remote_size_value" && "$remote_size_value" == "$local_size" ]]; then
      break
    fi
    if [[ "$attempt" -lt 4 ]]; then
      sleep 20
    fi
    attempt=$((attempt + 1))
  done

  if [[ "$remote_size_value" == "$local_size" ]]; then
    log_message "Uploaded {.file $rel}" --message-type success
  elif [[ -z "$remote_size_value" ]]; then
    log_message "Uploaded {.file $rel} (remote size unavailable; size check skipped)" --message-type warning
  else
    log_message "Size mismatch after upload: {.file $rel} remote={.val {$remote_size_value}} local={.val {$local_size}}" --message-type error
    exit 1
  fi
}

download_remote_md5sum() {
  local out="$1"
  local attempt=1
  while [[ "$attempt" -le 20 ]]; do
    if curl --silent --show-error --fail \
      --connect-timeout 60 \
      --max-time 300 \
      --user "${FTP_USER}:${FTP_PASSWORD}" \
      "ftp://${FTP_HOST}:${FTP_PORT}/${REMOTE_DIR}/md5sum.txt" \
      --output "$out"; then
      return 0
    fi
    log_message "Retrying remote {.file md5sum.txt} download (attempt {.val $attempt}/20)" --message-type warning
    sleep 60
    attempt=$((attempt + 1))
  done
  return 1
}

changed_files_by_md5() {
  local remote_md5="$1"
  local local_md5="md5sum.txt"

  if [[ ! -f "$local_md5" ]]; then
    log_message "Local {.file md5sum.txt} does not exist under {.path $LOCAL_DIR}" --message-type error
    exit 1
  fi

  awk '
    NR == FNR {
      if (NF >= 2) {
        md5 = $1
        file = $2
        for (i = 3; i <= NF; i++) file = file " " $i
        remote[file] = md5
      }
      next
    }
    NF >= 2 {
      md5 = $1
      file = $2
      for (i = 3; i <= NF; i++) file = file " " $i
      if (!(file in remote) || remote[file] != md5) print file
    }
  ' "$remote_md5" "$local_md5"
}

cd "$LOCAL_DIR"

FILES=()
if [[ "$CHANGED_BY_MD5" -eq 1 ]]; then
  REMOTE_MD5_FILE="$(mktemp)"
  trap 'rm -f "$REMOTE_MD5_FILE"' EXIT
  if download_remote_md5sum "$REMOTE_MD5_FILE"; then
    while IFS= read -r rel; do
      [[ -n "$rel" ]] && FILES+=("$rel")
    done < <(changed_files_by_md5 "$REMOTE_MD5_FILE" | grep -v '^md5sum\.txt$' || true)
    if [[ "${#FILES[@]}" -gt 0 ]] || changed_files_by_md5 "$REMOTE_MD5_FILE" | grep -qx 'md5sum.txt'; then
      FILES+=("md5sum.txt")
    fi
  else
    log_message "Could not download remote {.file md5sum.txt}; refusing md5-based upload to avoid unintended large transfers." --message-type error
    exit 1
  fi
fi

if [[ "${#FILES[@]}" -eq 0 && "$CHANGED_BY_MD5" -eq 0 ]]; then
  while IFS= read -r rel; do
    FILES+=("$rel")
  done < <(
    find . -type f |
      sed 's#^\./##' |
      awk '{ print ($0 == "expression/matrix.mtx.gz" ? "1\t" : "0\t") $0 }' |
      sort -k1,1n -k2,2 |
      cut -f2-
  )
fi

log_message "Uploading {.val {${#FILES[@]}}} files from {.path $LOCAL_DIR} to {.url ftp://${FTP_HOST}:${FTP_PORT}/${REMOTE_DIR}}" --message-type running

if [[ "${#FILES[@]}" -eq 0 ]]; then
  log_message "No files need upload." --message-type success
  exit 0
fi

if [[ "$DELETE_REMOTE_EXTRA" -eq 1 ]]; then
  log_message "Deleting known obsolete remote files under {.url ftp://${FTP_HOST}:${FTP_PORT}/${REMOTE_DIR}}" --message-type running
  for rel in \
    "embeddings/pca.tsv.gz" \
    "embeddings/umap.tsv.gz" \
    "validation/lisi_summary.tsv" \
    "provenance/reference_audit.tsv" \
    "provenance/source_access_summary.tsv" \
    "provenance/dataset_summary.tsv" \
    "figures/fig4_age_coverage.svg" \
    "figures/fig4_age_coverage.pdf" \
    "figures/fig1_revised.pdf"; do
    log_message "Deleting remote obsolete file {.file $rel}" --message-type warning
    delete_remote_file "$rel"
  done
  log_message "Removing remote empty directory {.file figures}" --message-type warning
  delete_remote_dir "figures"
fi

for rel in "${FILES[@]}"; do
  upload_file "$rel"
done

log_message "Upload complete." --message-type success
