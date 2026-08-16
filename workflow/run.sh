#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
MODE="${1:-verify}"
case "$MODE" in
  verify)
    exec python3 "$ROOT_DIR/tests/accept_release.py" "$ROOT_DIR"
    ;;
  render|analysis-frozen|full-raw)
    printf '%s\n' "$MODE: DOCUMENTED_NOT_YET_VALIDATED in Gate12BM" >&2
    exit 3
    ;;
  *)
    printf '%s\n' "usage: workflow/run.sh {verify|render|analysis-frozen|full-raw}" >&2
    exit 2
    ;;
esac
