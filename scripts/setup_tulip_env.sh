#!/usr/bin/env bash
# Build an isolated Python environment for running TULIP-TCR out-of-process,
# for use with mhctools' `Tulip` predictor (mhctools/tulip.py).
#
# Why isolated: TULIP-TCR is GPLv3 and pinned to transformers==4.32.1 (which in
# turn needs an old tokenizers). mhctools is Apache-2.0 and depends on neither
# torch nor transformers, so we never install those into your mhctools
# environment -- we build a *separate* interpreter and shell out to it.
#
# Why Python 3.11: transformers==4.32.1 resolves tokenizers 0.13.x, which has
# NO cp312 wheel (it predates Python 3.12) and would otherwise be built from
# source (requiring a Rust toolchain). Python 3.11 has a prebuilt wheel, so the
# install is wheels-only and needs no compiler.
#
# Usage:
#   scripts/setup_tulip_env.sh [ENV_DIR] [TULIP_HOME]
#
#   ENV_DIR     where to create the venv        (default: ./tulip-env)
#   TULIP_HOME  where to clone/find TULIP-TCR   (default: ./TULIP-TCR)
#
# After it runs, export the two variables it prints, e.g.:
#   export TULIP_HOME=/path/to/TULIP-TCR
#   export TULIP_PYTHON=/path/to/tulip-env/bin/python

set -euo pipefail

ENV_DIR="${1:-./tulip-env}"
TULIP_HOME="${2:-./TULIP-TCR}"
PY311="${TULIP_SETUP_PYTHON:-python3.11}"
TULIP_REPO="https://github.com/barthelemymp/TULIP-TCR.git"

log() { printf '[setup_tulip_env] %s\n' "$*" >&2; }

if ! command -v "$PY311" >/dev/null 2>&1; then
    log "ERROR: '$PY311' not found. Install Python 3.11 (brew install python@3.11)"
    log "or set TULIP_SETUP_PYTHON to a 3.11 interpreter."
    exit 1
fi

# Clone TULIP-TCR (GPLv3) if it isn't already present. We never redistribute
# it; you obtain it yourself under its own license.
if [[ ! -f "$TULIP_HOME/predict.py" ]]; then
    log "Cloning TULIP-TCR (GPLv3) into $TULIP_HOME"
    git clone --depth 1 "$TULIP_REPO" "$TULIP_HOME"
else
    log "Found existing TULIP-TCR at $TULIP_HOME"
fi

# Create the isolated env. Prefer uv (fast) but fall back to venv+pip.
DEPS=(torch "transformers==4.32.1" scikit-learn pandas numpy)
if command -v uv >/dev/null 2>&1; then
    log "Creating venv with uv at $ENV_DIR (Python 3.11)"
    uv venv --python "$PY311" "$ENV_DIR"
    log "Installing: ${DEPS[*]}"
    VIRTUAL_ENV="$ENV_DIR" uv pip install "${DEPS[@]}"
else
    log "uv not found; using venv + pip at $ENV_DIR"
    "$PY311" -m venv "$ENV_DIR"
    "$ENV_DIR/bin/python" -m pip install --upgrade pip
    log "Installing: ${DEPS[*]}"
    "$ENV_DIR/bin/python" -m pip install "${DEPS[@]}"
fi

ENV_PY="$(cd "$ENV_DIR" && pwd)/bin/python"
TULIP_ABS="$(cd "$TULIP_HOME" && pwd)"

# Smoke test: TULIP's model code must import in the isolated interpreter.
log "Verifying TULIP imports in the isolated interpreter ..."
( cd "$TULIP_ABS" && "$ENV_PY" -c "import src.multiTrans; import torch, transformers; \
print('ok: transformers', transformers.__version__, '| torch', torch.__version__)" )

cat <<EOF

Done. Add these to your shell / CI env to use mhctools.Tulip:

  export TULIP_HOME=$TULIP_ABS
  export TULIP_PYTHON=$ENV_PY

EOF
