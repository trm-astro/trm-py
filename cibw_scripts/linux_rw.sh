#!/bin/bash
set -euo pipefail
yum install -y jq

# Define the Conan package info for PCRE2
PCRE_ID=$(conan list 'pcre2/10.44:*' --format=json | jq -r '.[]."pcre2/10.44".revisions[].packages | keys[0]')
PCRE_LIB=$(conan cache path pcre2/10.44:$PCRE_ID)

# Ensure libpcre2-posix is correctly linked (manually or via RPATH/RUNPATH if needed)
LIB_FILE="$PCRE_LIB/lib/libpcre2-posix.so.3.0.5"

if ! ldd "$LIB_FILE" | grep -q "$PCRE_LIB/lib"; then
    echo "Setting RPATH to locate shared libraries"
    patchelf --set-rpath '$ORIGIN' "$LIB_FILE"
fi



#getting into the weeds
# find the wheel and extract the library
rm -rf uncomp || true
unzip "$2" -d "uncomp"
# inspect the three c++ extensions found in the wheel
# check the shared object files
echo "Checking trm_py/_cpp"
ls -l uncomp/trm_py/_cpp/

echo "Checking shared object file doppler"
auditwheel lddtree uncomp/trm_py/_cpp/_cpp_doppler*.so

echo "Checking shared object file roche"
auditwheel lddtree uncomp/trm_py/_cpp/_cpp_roche*.so

echo "Checking shared object file subs"
auditwheel lddtree uncomp/trm_py/_cpp/_cpp_subs*.so

# Repair the wheel using auditwheel
echo "Auditwheel show"
auditwheel show "$2"
echo "Auditwheel repair: auditwheel repair $2 -w $1"
auditwheel repair "$2" -w "$1"