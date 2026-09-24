#!/bin/bash

# Enhanced download utilities with resume capability, progress display, and integrity checking
# Usage: source this file in other download scripts

set -e

# Reuse the new unified shell logger.
UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$UTILS_DIR/log_message.sh"

file_size_bytes() {
    local file_path="$1"
    local size

    if size=$(stat -f%z "$file_path" 2>/dev/null); then
        printf '%s\n' "$size"
        return 0
    fi

    if size=$(stat -c%s "$file_path" 2>/dev/null); then
        printf '%s\n' "$size"
        return 0
    fi

    wc -c < "$file_path" | tr -d ' '
}

# Expose the DOWNLOAD_LIST declared by a dataset downloader as a read-only
# interface. This keeps URL, byte-size, and SHA-256 information in the same
# module that performs the download instead of duplicating it in another
# manifest.
validate_download_list() {
    local download_list="$1"
    local url filename expected_size expected_sha256 extra
    local row=0

    while IFS='|' read -r \
        url filename expected_size expected_sha256 extra; do
        [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
        row=$((row + 1))
        if [ -n "$extra" ] || [ -z "$filename" ]; then
            echo "Malformed download record at row $row" >&2
            return 1
        fi
        case "$filename" in
            */*|.|..|*[!A-Za-z0-9._-]*)
                echo "Unsafe download filename at row $row: $filename" >&2
                return 1
                ;;
        esac
        case "$expected_size" in
            ""|*[!0-9]*)
                echo "Invalid expected byte size at row $row: $expected_size" >&2
                return 1
                ;;
        esac
        if [ -n "$expected_sha256" ] &&
           ! [[ "$expected_sha256" =~ ^[0-9a-f]{64}$ ]]; then
            echo "Invalid SHA-256 at row $row for $filename" >&2
            return 1
        fi
    done <<< "$download_list"

    if [ "$row" -eq 0 ]; then
        echo "Download list contains no file records" >&2
        return 1
    fi
}

download_query() {
    local download_list="$1"
    shift
    local action="${1:-}"
    local requested_filename="${2:-}"
    local url filename expected_size expected_sha256
    local record=""
    local matches=0

    validate_download_list "$download_list"
    case "$action" in
        --list)
            if [ "$#" -ne 1 ]; then
                echo "Usage: <dataset downloader> --list" >&2
                return 2
            fi
            printf 'source_url|file_name|expected_bytes|sha256\n'
            while IFS='|' read -r \
                url filename expected_size expected_sha256; do
                [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
                printf '%s|%s|%s|%s\n' \
                    "$url" "$filename" "$expected_size" "$expected_sha256"
            done <<< "$download_list"
            ;;
        --record)
            if [ "$#" -ne 2 ] || [ -z "$requested_filename" ]; then
                echo \
                    "Usage: <dataset downloader> --record <filename>" \
                    >&2
                return 2
            fi
            while IFS='|' read -r \
                url filename expected_size expected_sha256; do
                [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
                if [ "$filename" = "$requested_filename" ]; then
                    record="$url|$filename|$expected_size|$expected_sha256"
                    matches=$((matches + 1))
                fi
            done <<< "$download_list"
            if [ "$matches" -ne 1 ]; then
                echo \
                    "Expected one download record for $requested_filename; found $matches" \
                    >&2
                return 1
            fi
            printf '%s\n' "$record"
            ;;
        *)
            echo \
                "Usage: <dataset downloader> --list | --record <filename>" \
                >&2
            return 2
            ;;
    esac
}

dataset_download_record() {
    local repo_dir="$1"
    local dataset="$2"
    local requested_filename="$3"
    local download_script="$repo_dir/download/$dataset.sh"
    local record url filename expected_size expected_sha256 extra

    if [ ! -f "$download_script" ]; then
        echo "Missing dataset downloader: $download_script" >&2
        return 1
    fi
    record="$(bash "$download_script" --record "$requested_filename")" ||
        return 1
    IFS='|' read -r \
        url filename expected_size expected_sha256 extra <<< "$record"
    if [ -n "$extra" ] || [ "$filename" != "$requested_filename" ]; then
        echo "Malformed download record for $dataset/$requested_filename" >&2
        return 1
    fi
    case "$expected_size" in
        ""|*[!0-9]*)
            echo \
                "Invalid expected byte size for $dataset/$requested_filename" \
                >&2
            return 1
            ;;
    esac
    if ! [[ "$expected_sha256" =~ ^[0-9a-f]{64}$ ]]; then
        echo "Missing or invalid SHA-256 for $dataset/$requested_filename" >&2
        return 1
    fi
    printf '%s|%s|%s|%s\n' \
        "$url" "$filename" "$expected_size" "$expected_sha256"
}

# Download implementation. Call it through download_with_resume(), which holds
# a per-output lock for the complete operation.
_download_with_resume_unlocked() {
    local url="$1"
    local output_file="$2"
    local max_retries="${3:-3}"
    local retry_delay="${4:-10}"
    local curl_max_time="${BRAINOMICS_CURL_MAX_TIME:-3600}"
    case "$curl_max_time" in
        ""|*[!0-9]*)
            log_message \
                "BRAINOMICS_CURL_MAX_TIME must be a non-negative integer" \
                --message-type error || true
            return 2
            ;;
    esac
    
    local filename=$(basename "$output_file")
    local temp_file="${output_file}.tmp"
    local resume_file="${output_file}.resume"
    
    # Create directory if it doesn't exist
    mkdir -p "$(dirname "$output_file")"
    
    # A resume marker means output_file is a saved partial download.
    if [ -f "$output_file" ] && [ ! -f "$resume_file" ]; then
        log_message "File {.file ${filename}} already exists, skipping download"
        return 0
    fi
    
    # Resume either an in-progress .tmp file or a partial output saved after a
    # failed attempt. curl -C - derives the offset from the file passed to -o.
    local resume_flag=""
    local resume_size=0
    if [ -f "$temp_file" ]; then
        resume_size=$(file_size_bytes "$temp_file")
    elif [ -f "$resume_file" ] && [ -f "$output_file" ]; then
        resume_size=$(file_size_bytes "$output_file")
    elif [ -f "$resume_file" ]; then
        rm -f "$resume_file"
    fi

    if [ "$resume_size" -gt 0 ]; then
        echo "$resume_size" > "$resume_file"
        resume_flag="-C -"
        log_message "Resuming download of {.file ${filename}} from byte {.val ${resume_size}}" --message-type warning
    elif [ -f "$output_file" ]; then
        local existing_size
        existing_size=$(file_size_bytes "$output_file")
        if [ "$existing_size" -gt 0 ]; then
            echo "$existing_size" > "$resume_file"
            resume_flag="-C -"
            log_message "Resuming download of {.file ${filename}} from byte {.val ${existing_size}}" --message-type warning
        fi
    fi
    
    local attempt=1
    while [ $attempt -le $max_retries ]; do
        log_message "Downloading {.file ${filename}} (attempt {.val ${attempt}}/{.val ${max_retries}})..."
        
        # Download with clean output and lightweight progress (every 10s)
        local attempt_status=0
        (
            # Obtain total size if available for percentage display
            total_size=$(curl -sI "$url" | awk -F": " 'tolower($1)=="content-length"{print $2}' | tr -d '\r')
            [ -z "$total_size" ] && total_size=0
            connecting_logged=0
            last_logged_mb=-1

            # If resuming, start from the resume position
            if [ -n "$resume_flag" ] && [ -f "$resume_file" ]; then
                resume_size=$(cat "$resume_file")
                if [ "$resume_size" -gt 0 ]; then
                    # Copy existing partial file to temp file for resume
                    if [ -f "$output_file" ] && [ ! -f "$temp_file" ]; then
                        cp "$output_file" "$temp_file"
                    fi
                fi
            fi

            curl -L \
                --connect-timeout 30 \
                --max-time "$curl_max_time" \
                --retry 3 \
                --retry-delay 10 \
                --retry-max-time 300 \
                --fail \
                --show-error \
                --silent \
                $resume_flag \
                -o "$temp_file" \
                "$url" &
            curl_pid=$!
            while kill -0 "$curl_pid" 2>/dev/null; do
                if [ -f "$temp_file" ]; then
                    current_size=$(file_size_bytes "$temp_file")
                    if [ "$current_size" -eq 0 ]; then
                        if [ "$connecting_logged" -eq 0 ]; then
                            log_message "Downloading {.file ${filename}} ... connecting" --message-type running
                            connecting_logged=1
                        fi
                    else
                        current_mb=$((current_size / 1024 / 1024))
                        if [ "$current_mb" -ne "$last_logged_mb" ]; then
                            if [ "$total_size" -gt 0 ]; then
                                percent=$(( current_size * 100 / total_size ))
                                log_message "Downloading {.file ${filename}} ... {.pkg ${current_mb}MB} ({.pkg ${percent}%})" --message-type running
                            else
                                log_message "Downloading {.file ${filename}} ... {.pkg ${current_mb}MB}" --message-type running
                            fi
                            last_logged_mb=$current_mb
                        fi
                    fi
                fi
                sleep 10
            done
            wait "$curl_pid"
        ) || attempt_status=$?
        if [ "$attempt_status" -eq 0 ]; then
            # Download successful
            mv "$temp_file" "$output_file"
            rm -f "$resume_file"
            
            # Display file size information
            local downloaded_size
            downloaded_size=$(file_size_bytes "$output_file")
            local size_mb=$((downloaded_size / 1024 / 1024))
            log_message "Successfully downloaded {.file ${filename}} ({.val ${size_mb}MB})" --message-type success
            return 0
        else
            local exit_code="$attempt_status"
            log_message "Download failed for {.file ${filename}} (exit code: {.val ${exit_code}})" --message-type error || true
            
            if [ -f "$temp_file" ]; then
                local current_size
                current_size=$(file_size_bytes "$temp_file")
                if [ "$current_size" -gt 0 ]; then
                    mv "$temp_file" "$output_file"
                    echo "$current_size" > "$resume_file"
                    resume_flag="-C -"
                    log_message "Saved resume position: {.val ${current_size}} bytes"
                else
                    rm -f "$temp_file"
                fi
            fi
            
            if [ $attempt -lt $max_retries ]; then
                log_message "Waiting {.val ${retry_delay}} seconds before retry..." --message-type warning
                sleep $retry_delay
            fi
        fi
        
        attempt=$((attempt + 1))
    done
    
    log_message "Failed to download {.file ${filename}} after {.pkg ${max_retries}} attempts" --message-type error || true
    rm -f "$temp_file"
    return 1
}

# Download function with resume capability, progress display, and an atomic
# per-output lock. The lock prevents two hosts or shell sessions that share a
# filesystem from writing the same target and temporary file concurrently.
download_with_resume() (
    local output_file="$2"
    local lock_dir="${output_file}.download.lock"

    mkdir -p "$(dirname "$output_file")"
    if ! mkdir "$lock_dir" 2>/dev/null; then
        log_message "Another download already owns {.file $(basename "$output_file")}; lock: {.file ${lock_dir}}" --message-type error || true
        return 75
    fi
    cleanup_download_lock() {
        rmdir "$lock_dir" 2>/dev/null || true
    }
    trap cleanup_download_lock EXIT
    trap 'exit 130' HUP INT TERM

    local status=0
    _download_with_resume_unlocked "$@" || status=$?
    return "$status"
)

# Verify file integrity using file size and basic checks
verify_file_integrity() {
    local file_path="$1"
    local expected_size="${2:-0}"
    local expected_sha256="${3:-}"
    
    if [ ! -f "$file_path" ]; then
        log_message "File $file_path does not exist" --message-type error || true
        return 1
    fi
    
    local actual_size
    actual_size=$(file_size_bytes "$file_path")
    
    if [ "$expected_size" -gt 0 ] && [ "$actual_size" -ne "$expected_size" ]; then
        log_message "File $file_path has the wrong byte size (expected: $expected_size, actual: $actual_size)" --message-type error || true
        return 1
    fi

    if [ -n "$expected_sha256" ]; then
        local actual_sha256
        if command -v sha256sum >/dev/null 2>&1; then
            actual_sha256=$(sha256sum "$file_path" | awk '{print $1}')
        elif command -v shasum >/dev/null 2>&1; then
            actual_sha256=$(shasum -a 256 "$file_path" | awk '{print $1}')
        else
            log_message \
                "No SHA-256 command is available for $file_path" \
                --message-type error || true
            return 1
        fi
        if [ "$actual_sha256" != "$expected_sha256" ]; then
            log_message \
                "File $file_path has the wrong SHA-256 checksum" \
                --message-type error || true
            return 1
        fi
    fi
    
    # Basic file type check for common formats
    local filename=$(basename "$file_path")
    case "$filename" in
        *.gz|*.bgz)
            if ! gzip -t "$file_path" 2>/dev/null; then
                log_message "File $file_path is not a valid gzip/BGZF file" --message-type error || true
                return 1
            fi
            ;;
        *.h5ad)
            local python_bin="${BRAINOMICS_PYTHON:-python3}"
            if command -v "$python_bin" >/dev/null 2>&1 &&
               "$python_bin" -c "import h5py" >/dev/null 2>&1; then
                if ! "$python_bin" - "$file_path" <<'PY'
import sys

import h5py

with h5py.File(sys.argv[1], "r") as handle:
    required = {"X", "obs", "var"}
    missing = required.difference(handle.keys())
    if missing:
        raise SystemExit(
            "missing required H5AD groups: " + ", ".join(sorted(missing))
        )
PY
                then
                    log_message "File $file_path is not a readable H5AD file" --message-type error || true
                    return 1
                fi
            else
                local hdf5_magic
                hdf5_magic=$(od -An -tx1 -N8 "$file_path" | tr -d ' \n')
                if [ "$hdf5_magic" != "894844460d0a1a0a" ]; then
                    log_message "File $file_path does not have an HDF5 signature" --message-type error || true
                    return 1
                fi
            fi
            ;;
        *.xlsx)
            if command -v unzip >/dev/null 2>&1 &&
               ! unzip -tqq "$file_path" >/dev/null 2>&1; then
                log_message "File $file_path is not a valid XLSX archive" --message-type error || true
                return 1
            fi
            ;;
    esac
    
    log_message "File integrity check passed for $filename" --message-type success
    return 0
}

batch_download() {
    local download_list="$1"
    local data_dir="$2"
    local max_retries="${3:-3}"
    
    log_message "Starting batch download of $(echo "$download_list" | wc -l) files..."
    
    local success_count=0
    local total_count=0
    
    while IFS='|' read -r url filename expected_size expected_sha256; do
        # Skip empty lines and comments
        [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
        
        total_count=$((total_count + 1))
        local output_file="$data_dir/$filename"
        
        log_message "Processing file {.val ${total_count}}: {.file ${filename}}"
        
        # A resume marker proves this is an intentionally incomplete file.
        # Resume it directly instead of spending minutes validating a large
        # gzip stream that cannot yet pass the complete-file gate.
        if [ -f "$output_file" ]; then
            if [ -f "${output_file}.resume" ]; then
                log_message "Preserving partial {.file ${filename}} for resume" --message-type warning
            elif verify_file_integrity \
                "$output_file" "$expected_size" "$expected_sha256"; then
                success_count=$((success_count + 1))
                echo ""
                continue
            else
                log_message "Integrity check failed for {.file ${filename}}, re-downloading..." --message-type error || true
                rm -f "$output_file"
            fi
        fi

        if download_with_resume "$url" "$output_file" "$max_retries"; then
            if verify_file_integrity \
                "$output_file" "$expected_size" "$expected_sha256"; then
                success_count=$((success_count + 1))
            else
                log_message "Integrity check failed for {.file ${filename}}" --message-type error || true
            fi
        else
            log_message "Download failed for {.file ${filename}}" --message-type error || true
        fi
        
        echo "" # Add blank line for readability
    done <<< "$download_list"
    
    log_message "Batch download completed: {.val ${success_count}}/{.val ${total_count}} files successful"
    
    if [ $success_count -eq $total_count ]; then
        log_message "All files downloaded successfully!" --message-type success
        return 0
    else
        log_message "Some files failed to download" --message-type error || true
        return 1
    fi
}

# Parallel batch download with concurrency control
batch_download_parallel() {
    local download_list="$1"
    local data_dir="$2"
    local max_retries="${3:-3}"
    local concurrency="${4:-3}"

    log_message "Starting parallel batch download (concurrency={.val ${concurrency}}) of {.val $(echo "$download_list" | wc -l)} files..."

    local temp_dir="$data_dir/.dl_tmp_$$"
    mkdir -p "$temp_dir"

    # Helper to process a single file (runs in background)
    _perform_single_download() {
        local url="$1"
        local filename="$2"
        local expected_size="$3"
        local expected_sha256="$4"
        local data_dir="$5"
        local max_retries="$6"

        local output_file="$data_dir/$filename"
        local status_file="$temp_dir/${filename}.status"

        # Skip complete-file validation for an explicitly resumable partial.
        if [ -f "$output_file" ]; then
            if [ -f "${output_file}.resume" ]; then
                log_message "Preserving partial {.file ${filename}} for resume" --message-type warning
            elif verify_file_integrity \
                "$output_file" "$expected_size" "$expected_sha256"; then
                echo ok > "$status_file"
                return 0
            else
                rm -f "$output_file"
            fi
        fi

        if download_with_resume "$url" "$output_file" "$max_retries"; then
            if verify_file_integrity \
                "$output_file" "$expected_size" "$expected_sha256"; then
                echo ok > "$status_file"
            else
                echo fail > "$status_file"
            fi
        else
            echo fail > "$status_file"
        fi
    }

    # Launch jobs with concurrency control
    local -a pids=()
    local total_count=0

    while IFS='|' read -r url filename expected_size expected_sha256; do
        [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
        total_count=$((total_count + 1))

        # Throttle to concurrency
        while [ ${#pids[@]} -ge $concurrency ]; do
            wait "${pids[0]}" 2>/dev/null || true
            pids=("${pids[@]:1}")
        done

        log_message "Queueing file {.val ${total_count}}: {.file ${filename}}"
        _perform_single_download \
            "$url" "$filename" "$expected_size" "$expected_sha256" \
            "$data_dir" "$max_retries" &
        pids+=("$!")
    done <<< "$download_list"

    # Wait for remaining jobs
    for pid in "${pids[@]}"; do
        wait "$pid" 2>/dev/null || true
    done

    # Summarize results
    local success_count=$(grep -l "^ok$" "$temp_dir"/*.status 2>/dev/null | wc -l | tr -d ' ')
    local total_status=$(ls "$temp_dir"/*.status 2>/dev/null | wc -l | tr -d ' ')
    [ -z "$total_status" ] && total_status=0

    log_message "Parallel batch download completed: {.val ${success_count}}/{.val ${total_status}} files successful"

    rm -rf "$temp_dir"

    if [ "$success_count" -eq "$total_status" ]; then
        log_message "All files downloaded successfully!" --message-type success
        return 0
    else
        log_message "Some files failed to download" --message-type error || true
        return 1
    fi
}


cleanup_temp_files() {
    local data_dir="$1"
    log_message "Cleaning up temporary files..."
    find "$data_dir" -name "*.tmp" -delete 2>/dev/null || true
    find "$data_dir" -name "*.resume" -delete 2>/dev/null || true
    log_message "Cleanup completed" --message-type success
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_message "$1 is not installed. Please install it first." --message-type error || true
        exit 1
    fi
}

run_r_script() {
    local base_dir="$1"
    local script_name="$2"
    local description="$3"
    shift 3
    local script_path="${base_dir}/${script_name}"
    if [ ! -f "$script_path" ]; then
        log_message "Script not found: $script_path" --message-type error || true
        exit 1
    fi
    if "${BRAINOMICS_RSCRIPT:-Rscript}" "$script_path" "$@"; then
        log_message "$description completed" --message-type success
    else
        log_message "$description failed" --message-type error || true
        exit 1
    fi
    echo ""
}

run_python_script() {
    local base_dir="$1"
    local script_name="$2"
    local description="$3"
    local run_from_script_dir="${4:-}"
    local script_path="${base_dir}/${script_name}"
    if [ ! -f "$script_path" ]; then
        log_message "Script not found: $script_path" --message-type error || true
        exit 1
    fi
    if [ -n "$run_from_script_dir" ]; then
        if (cd "$(dirname "$script_path")" && "${BRAINOMICS_PYTHON:-python3}" "$(basename "$script_path")"); then
            log_message "$description completed" --message-type success
        else
            log_message "$description failed" --message-type error || true
            exit 1
        fi
    else
        if "${BRAINOMICS_PYTHON:-python3}" "$script_path"; then
            log_message "$description completed" --message-type success
        else
            log_message "$description failed" --message-type error || true
            exit 1
        fi
    fi
    echo ""
}
