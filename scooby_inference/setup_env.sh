#!/usr/bin/env bash
# Step 2a: build the inference environment on pasteur (run once, ~10 min).
# Uses uv (fast pip replacement). Creates ~/venvs/scooby-inf with Python 3.11.
set -euo pipefail
cd "$(dirname "$0")"

command -v uv >/dev/null || python3 -m pip install --user uv
export PATH="$HOME/.local/bin:$PATH"

uv venv ~/venvs/scooby-inf -p 3.11
source ~/venvs/scooby-inf/bin/activate

# PyTorch built for CUDA 12.8: the first wheels with Blackwell (sm_120) kernels.
uv pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128

# Everything else. torch is left out on purpose so the CUDA 12.8 build above is kept.
uv pip install -r requirements-step2.txt

# scooby itself without its declared dependencies: its pins (pybigtools==0.1.1, open-ended
# borzoi/enformer versions) do not resolve together on Python 3.11, and pybigtools is only
# needed for training.
uv pip install --no-deps git+https://github.com/gagneurlab/scooby.git

# The authors' variant-scoring script and gene windows.
[ -d scooby_reproducibility ] || git clone --depth 1 https://github.com/gagneurlab/scooby_reproducibility.git

echo "Environment ready. Next: python smoke_test.py"
