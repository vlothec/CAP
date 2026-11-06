#!/bin/bash
# setup_conda.sh
# Description: Creates the conda environment and installs the BCT package.

set -e

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ENV_NAME="${CAP_ENV_NAME:-cap-pipeline}"
ENV_FILE="environment.yml"
LOCK_FILE="conda-linux-64.lock"

# Create the environment first so setup does not require Git on the host.
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "Conda environment '$ENV_NAME' already exists; reusing it."
elif [ "$(uname -s)" = "Linux" ] && [ "$(uname -m)" = "x86_64" ] && [ -f "$LOCK_FILE" ]; then
    echo "Creating Conda environment '$ENV_NAME' from reproducible Linux lock..."
    conda create --name "$ENV_NAME" --file "$LOCK_FILE"
else
    echo "No matching lock for this platform; creating '$ENV_NAME' from $ENV_FILE..."
    conda env create --name "$ENV_NAME" -f "$ENV_FILE"
fi

# Use Git installed inside the Conda environment to initialize TRASH2.
if [ -z "$(ls -A modules/TRASH_2 2>/dev/null)" ]; then
    echo "Initializing submodules..."
    conda run -n "$ENV_NAME" git submodule update --init --recursive
fi

echo "Activating environment..."
# Need to source conda.sh to use 'conda activate' in script, or use 'conda run'
# Assuming 'conda run' is available (newer conda versions)

echo "Setting up BCT (Bayesian Context Trees)..."
echo "Attempting to compile local C++ binary (faster)..."

# Try compilation and an execution smoke test first
if CAP_ENV_NAME="$ENV_NAME" CAP_CONDA_EXE=conda scripts/build-ctw.sh; then
    echo "✓ C++ binary compiled and executed successfully."
else
    echo "⚠️  C++ build or execution validation failed. Falling back to CRAN installation (slower)..."
    conda run -n "$ENV_NAME" Rscript install_bioc_packages.R
fi

echo "Setup complete! Activate the environment with:"
echo "conda activate $ENV_NAME"
