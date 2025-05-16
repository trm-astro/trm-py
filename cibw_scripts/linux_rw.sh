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

# Repair the wheel using auditwheel
auditwheel show "$2"
auditwheel repair "$2" -w "$1"