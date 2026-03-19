#!/bin/bash
set -euo pipefail

# Compile Cython extensions in-place
cd "${SRC_DIR}/src"
python setup.py build_ext --inplace

# Install .so files directly into site-packages so bare 'import calc' etc. resolves from anywhere
cp calc*.so parse*.so thread_wrapper*.so "${SP_DIR}/"

# Install scripts to bin/ so they are on PATH
cp "${SRC_DIR}/src/afcn.py" "${PREFIX}/bin/afcn.py"
cp "${SRC_DIR}/src/convert_genexpc_to_afcn.py" "${PREFIX}/bin/convert_genexpc_to_afcn.py"
chmod +x "${PREFIX}/bin/afcn.py" "${PREFIX}/bin/convert_genexpc_to_afcn.py"
