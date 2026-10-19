#!/bin/zsh
# Double-click in Finder to open the local app; no notebook kernel is required.
set -eu
PROJECT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
cd "$PROJECT_DIR"
APP_FILE="$PROJECT_DIR/outputs/spatial_handoff/app/index.html"
if [[ -x "$PROJECT_DIR/.venv/bin/python" ]]; then
  PYTHON_BIN="$PROJECT_DIR/.venv/bin/python"
else
  PYTHON_BIN="python3"
fi
# Some Macs have a python3 stub but no usable Python installation. The saved
# app can still open directly, so do not make a Python install a viewing step.
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1 || ! "$PYTHON_BIN" -S -c 'import sys; raise SystemExit(sys.version_info.major != 3)' >/dev/null 2>&1; then
  if [[ -f "$APP_FILE" ]]; then
    printf 'Python is unavailable; opening the portable spatial app directly.\n'
    printf 'Notebook localhost links need a running server; this file has the same saved results.\n'
    open "$APP_FILE"
    exit 0
  fi
  printf 'The saved app is missing: %s\nDownload the complete repository to view saved results.\n' "$APP_FILE"
  exit 1
fi
if ! "$PYTHON_BIN" "$PROJECT_DIR/scripts/serve_spatial_app.py" "$@"; then
  printf '\nPress Return to close this window.\n'
  read -r REPLY
  exit 1
fi
