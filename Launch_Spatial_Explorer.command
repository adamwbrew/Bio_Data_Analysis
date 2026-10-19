#!/bin/zsh
# Double-click in Finder to open the local app; no notebook kernel is required.
set -eu
PROJECT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
cd "$PROJECT_DIR"
if [[ -x "$PROJECT_DIR/.venv/bin/python" ]]; then
  PYTHON_BIN="$PROJECT_DIR/.venv/bin/python"
else
  PYTHON_BIN="python3"
fi
if ! "$PYTHON_BIN" "$PROJECT_DIR/scripts/serve_spatial_app.py" "$@"; then
  printf '\nPress Return to close this window.\n'
  read -r REPLY
  exit 1
fi
