#!/usr/bin/env bash

# Shared resume/checkpoint helpers for SEN smoke runners.
# A stage is skipped only when:
#   1) its validator confirms all expected outputs are complete, and
#   2) its checkpoint signature matches the current scripts/config.
# Existing validated outputs without a checkpoint can be adopted once so that
# introducing checkpointing does not force needless recomputation.

CHECKPOINT_ROOT="${SEN_CHECKPOINT_ROOT:-${SEN_ROOT:?SEN_ROOT must be set}/.checkpoints}"
mkdir -p "$CHECKPOINT_ROOT"

checkpoint_signature() {
  local f
  {
    for f in "$@"; do
      if [[ -f "$f" ]]; then
        sha256sum "$f"
      else
        printf 'MISSING  %s\n' "$f"
      fi
    done
  } | sha256sum | awk '{print $1}'
}

checkpoint_should_skip() {
  local stage="$1"
  local validator="$2"
  shift 2
  local sig stamp old

  if [[ "${SEN_FORCE_RERUN:-0}" == "1" ]]; then
    echo "[RUN]  $stage: forced by SEN_FORCE_RERUN=1"
    return 1
  fi

  if ! "$validator"; then
    echo "[RUN]  $stage: expected outputs are incomplete or invalid"
    return 1
  fi

  sig="$(checkpoint_signature "$@")"
  stamp="$CHECKPOINT_ROOT/$stage.sha256"

  if [[ -s "$stamp" ]]; then
    old="$(tr -d '[:space:]' < "$stamp")"
    if [[ "$old" == "$sig" ]]; then
      echo "[SKIP] $stage: validated outputs + matching checkpoint"
      return 0
    fi
    echo "[RUN]  $stage: workflow/config signature changed"
    return 1
  fi

  if [[ "${SEN_ADOPT_EXISTING_OUTPUTS:-1}" == "1" ]]; then
    printf '%s\n' "$sig" > "$stamp"
    echo "[SKIP] $stage: validated pre-existing outputs; checkpoint adopted"
    return 0
  fi

  echo "[RUN]  $stage: outputs exist but no checkpoint is registered"
  return 1
}

checkpoint_mark() {
  local stage="$1"
  shift
  checkpoint_signature "$@" > "$CHECKPOINT_ROOT/$stage.sha256"
  echo "[CHECKPOINT] $stage"
}

checkpoint_require_valid() {
  local stage="$1"
  local validator="$2"
  if ! "$validator"; then
    echo "[FAIL] $stage completed but output validation failed." >&2
    exit 1
  fi
}
