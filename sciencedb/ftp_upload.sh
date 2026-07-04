#!/usr/bin/env bash
set -uo pipefail

BASE=${SCIDB_UPLOAD_BASE:-/Users/mx/Study/data/BrainOmicsData/ScienceDB}
REMOTE_DIR=${SCIDB_UPLOAD_REMOTE_DIR:-human_brain_age_interval_sc_snRNAseq_dataset}
HOST=${SCIDB_UPLOAD_HOST:-ftp-upload.scidb.cn}
PORT=${SCIDB_UPLOAD_PORT:-2121}
LOG=${SCIDB_UPLOAD_LOG:-/Users/mx/Study/data/BrainOmicsData/sciencedb_ftp_upload.log}
NETRC=${SCIDB_FTP_NETRC:?SCIDB_FTP_NETRC is required}

cd "$BASE" || exit 1

{
  echo "[$(date '+%F %T')] FTP upload started"
  echo "Base: $BASE"
  echo "Remote: ftp://$HOST:$PORT/$REMOTE_DIR/"
} >> "$LOG"

curl --silent --show-error --connect-timeout 30 --max-time 120 \
  --netrc-file "$NETRC" \
  -Q "DELE /$REMOTE_DIR/README.retry.md" \
  "ftp://$HOST:$PORT/" >> "$LOG" 2>&1 || true

find . -type f ! -name ".DS_Store" | sed "s#^./##" | sort | while IFS= read -r rel; do
  size=$(stat -f "%z" "$rel")
  echo "[$(date '+%F %T')] START $rel $size bytes" >> "$LOG"
  curl --silent --show-error --fail \
    --connect-timeout 60 \
    --retry 5 \
    --retry-delay 30 \
    --retry-all-errors \
    --ftp-create-dirs \
    --netrc-file "$NETRC" \
    -T "$rel" \
    "ftp://$HOST:$PORT/$REMOTE_DIR/$rel" >> "$LOG" 2>&1
  code=$?
  if [ "$code" -ne 0 ]; then
    echo "[$(date '+%F %T')] ERROR $rel curl_exit=$code" >> "$LOG"
    exit "$code"
  fi
  echo "[$(date '+%F %T')] DONE  $rel" >> "$LOG"
done

echo "[$(date '+%F %T')] FTP upload completed" >> "$LOG"
