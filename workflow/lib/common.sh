#!/usr/bin/env bash
# Shared path/config loader for SEN_genomics_pipeline.
# Source this near the top of shell workflow scripts.
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[1]:-${BASH_SOURCE[0]}}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

PATHS_FILE="${SEN_PATHS_FILE:-$REPO_ROOT/config/paths.env}"
if [[ -f "$PATHS_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$PATHS_FILE"
fi

export SEN_ROOT="${SEN_ROOT:-$REPO_ROOT}"
export SEN_RAW_READS="${SEN_RAW_READS:-$SEN_ROOT/fastq}"
export SEN_TRIMMED_READS="${SEN_TRIMMED_READS:-$SEN_ROOT/trimmed_fastq}"
export SEN_RESULTS="${SEN_RESULTS:-$SEN_ROOT/results}"
export SEN_REFERENCE_DIR="${SEN_REFERENCE_DIR:-$SEN_ROOT/reference}"

require_file() {
  [[ -s "$1" ]] || { echo "[ERROR] Missing/empty file: $1" >&2; exit 1; }
}

require_dir() {
  [[ -d "$1" ]] || { echo "[ERROR] Missing directory: $1" >&2; exit 1; }
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Required command not found: $1" >&2; exit 1; }
}
