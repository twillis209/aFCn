#!/bin/bash
set -euo pipefail

# Compile Cython extensions in-place
cd "${SRC_DIR}/src"
python setup.py build_ext --inplace

# Install afcn.py and compiled .so files as a package in site-packages
AFCN_PKG="${SP_DIR}/afcn"
mkdir -p "${AFCN_PKG}"

cp afcn.py calc*.so parse*.so thread_wrapper*.so "${AFCN_PKG}/"

# __init__.py so Python treats it as a package (not strictly needed, but tidy)
touch "${AFCN_PKG}/__init__.py"

# Entry-point script: cd into the package dir so bare 'import calc' resolves
mkdir -p "${PREFIX}/bin"
cat > "${PREFIX}/bin/afcn" << 'EOF'
#!/bin/bash
exec python -c "
import sys, os
pkg = os.path.join(sys.prefix, 'lib', 'python' + '.'.join(map(str, sys.version_info[:2])), 'site-packages', 'afcn')
sys.path.insert(0, pkg)
os.chdir(pkg)
exec(open(os.path.join(pkg, 'afcn.py')).read())
" -- "$@"
EOF
chmod +x "${PREFIX}/bin/afcn"
