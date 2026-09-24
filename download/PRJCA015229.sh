#!/bin/bash

# Download script for PRJCA015229 public NGDC OMIX files.
# BioProject: https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA015229
# OMIX005925: https://ngdc.cncb.ac.cn/omix/release/OMIX005925
# OMIX005926: https://ngdc.cncb.ac.cn/omix/release/OMIX005926
# OMIX005931: https://ngdc.cncb.ac.cn/omix/release/OMIX005931
# OMIX005932: https://ngdc.cncb.ac.cn/omix/release/OMIX005932

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${BRAINOMICS_DATA_ROOT:-$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData}"
DATA_DIR="$DATA_ROOT/raw/PRJCA015229"

log_message "Starting PRJCA015229 data download..."

DOWNLOAD_LIST="
ftp://download.cncb.ac.cn/OMIX/OMIX005925/OMIX005925-01.tar.gz|OMIX005925/OMIX005925-01.tar.gz|331916624
ftp://download.cncb.ac.cn/OMIX/OMIX005925/OMIX005925-02.tar.gz|OMIX005925/OMIX005925-02.tar.gz|257814162
ftp://download.cncb.ac.cn/OMIX/OMIX005925/OMIX005925-03.tar.gz|OMIX005925/OMIX005925-03.tar.gz|268732662
ftp://download.cncb.ac.cn/OMIX/OMIX005925/OMIX005925-04.tar.gz|OMIX005925/OMIX005925-04.tar.gz|349263904
ftp://download.cncb.ac.cn/OMIX/OMIX005925/OMIX005925-05.tar.gz|OMIX005925/OMIX005925-05.tar.gz|15197652
ftp://download.cncb.ac.cn/OMIX/OMIX005925/OMIX005925-06.tar.gz|OMIX005925/OMIX005925-06.tar.gz|67981049
ftp://download.cncb.ac.cn/OMIX/OMIX005926/OMIX005926-01.tar.gz|OMIX005926/OMIX005926-01.tar.gz|469464621
ftp://download.cncb.ac.cn/OMIX/OMIX005926/OMIX005926-02.tar.gz|OMIX005926/OMIX005926-02.tar.gz|36211203
ftp://download.cncb.ac.cn/OMIX/OMIX005926/OMIX005926-03.tar.gz|OMIX005926/OMIX005926-03.tar.gz|179646510
ftp://download.cncb.ac.cn/OMIX/OMIX005926/OMIX005926-04.tar.gz|OMIX005926/OMIX005926-04.tar.gz|259070987
ftp://download.cncb.ac.cn/OMIX/OMIX005926/OMIX005926-05.tar.gz|OMIX005926/OMIX005926-05.tar.gz|43429262
ftp://download.cncb.ac.cn/OMIX/OMIX005931/OMIX005931-01.tar|OMIX005931/OMIX005931-01.tar|1161463808
ftp://download.cncb.ac.cn/OMIX/OMIX005932/OMIX005932-01.txt|OMIX005932/OMIX005932-01.txt|38199625
ftp://download.cncb.ac.cn/OMIX/OMIX005932/OMIX005932-02.tar|OMIX005932/OMIX005932-02.tar|1005137920
"

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

for archive in "$DATA_DIR"/OMIX*/*.tar "$DATA_DIR"/OMIX*/*.tar.gz "$DATA_DIR"/OMIX*/*.tgz; do
    [ -f "$archive" ] || continue
    marker="$DATA_DIR/.${archive#"$DATA_DIR"/}_extracted"
    marker="${marker//\//_}"
    if [ ! -f "$marker" ]; then
        log_message "Extracting {.file ${archive#"$DATA_DIR"/}}..."
        tar -xf "$archive" -C "$DATA_DIR"
        touch "$marker"
    else
        log_message "Archive already extracted: {.file ${archive#"$DATA_DIR"/}}"
    fi
done

cleanup_temp_files "$DATA_DIR"

log_message "PRJCA015229 data download completed!" --message-type success
