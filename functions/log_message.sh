#!/usr/bin/env bash
# Shell entry point that reuses the log_message.sh shipped with thisutils.
#
#   source functions/log_message.sh
#   log_message "message" --message-type info
#
# Resolution order:
#   1. $LOG_MESSAGE_SH (explicit path to thisutils/log_message.sh)
#   2. R lookup: Rscript -e 'cat(system.file("scripts/log_message.sh", package="thisutils"))'
#   3. Plain fallback with the same call signature and ignored options
#
# Progress and information messages only; stdout parsed by other scripts (for
# example key<TAB>value lines) is left untouched.

_thisutils_log_message_source() {
  if [ -n "${LOG_MESSAGE_SH:-}" ] && [ -f "${LOG_MESSAGE_SH}" ]; then
    printf '%s\n' "${LOG_MESSAGE_SH}"; return 0
  fi
  if command -v "${BRAINOMICS_RSCRIPT:-Rscript}" >/dev/null 2>&1; then
    local resolved
    # Current layout: thisutils/scripts/; earlier installs used python/.
    resolved=$("${BRAINOMICS_RSCRIPT:-Rscript}" --vanilla -e 'p <- system.file("scripts/log_message.sh", package = "thisutils"); if (!nzchar(p)) p <- system.file("python/log_message.sh", package = "thisutils"); cat(p)' 2>/dev/null)
    if [ -n "$resolved" ] && [ -f "$resolved" ]; then printf '%s\n' "$resolved"; return 0; fi
  fi
  return 1
}

if _src=$(_thisutils_log_message_source); then
  # shellcheck source=/dev/null
  source "$_src"
else
  log_message() {
    local args=()
    while [ $# -gt 0 ]; do
      case "$1" in
        --message-type|--timestamp|--verbose|--cli-model|--indent|--color|--bg-color)
          shift 2 || shift ;;
        *) args+=("$1"); shift ;;
      esac
    done
    printf '%s\n' "${args[*]}" >&2
  }
  printf '%s\n' "log_message: thisutils not found, falling back to plain output (set LOG_MESSAGE_SH)" >&2
fi
