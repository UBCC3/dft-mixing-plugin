#!/bin/bash

echo "========= DFT Mixing Environment Setup ========"


echo "========= Installing Prerequisites ========"
conda install -c conda-forge \
    sqlalchemy \
    sqlite3 \
    pylibxc \
    psi4 \

PSIPATH="/usr/local/psi4/lib"

# In case user runs this from else where
cd "$(dirname "$0")"

# Check PSI4 location
PSIPATH=$(realpath "$(which psi4)/../../lib")
PLUGINPATH="${PWD}"

PSI_AND_MODULE_PATH=${PSIPATH}:${PWD}

export PYTHONPATH=${PSI_AND_MODULE_PATH}:${PYTHONPATH}
echo "Do you want to make a permanent change to your conda environment? [Y|n]"
read -p "> " user_resp

if [[ "$user_resp" =~ ^[Yy]$ || -z "$user_resp" ]]; then
    if [[ -z "$CONDA_PREFIX" ]]; then
        echo "No conda environment is active. Please 'conda activate <env>' before running this."
        cd "$prev_path"
        exit 1
    fi

    # Create Conda activation script
    mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
    cat <<EOF > "$CONDA_PREFIX/etc/conda/activate.d/dft_mixing_env.sh"
# Added by DFT Mixing installer
export PYTHONPATH="\$PYTHONPATH:$PSI_AND_MODULE_PATH"
EOF

    echo "Installed PYTHONPATH to Conda env: $CONDA_PREFIX"
else
    echo "Skipped making changes permanent."
fi

echo "Installation finished."

