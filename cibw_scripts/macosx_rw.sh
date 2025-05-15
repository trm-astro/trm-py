#! /bin/bash

#! define the cibw repair wheel script for macosx
# Relink the broken conan packages
PCRE_ID=$(conan list 'pcre2/10.44:*' --format=json | jq -r '.[]."pcre2/10.44".revisions[].packages | keys[0]')
PCRE_LIB=$(conan cache path pcre2/10.44:$PCRE_ID)
if ! otool -l $PCRE_LIB/lib/libpcre2-posix.3.0.5.dylib | grep -A2 LC_RPATH | grep -q "@loader_path"; then
    install_name_tool -add_rpath @loader_path $PCRE_LIB/lib/libpcre2-posix.3.0.5.dylib
fi

# Repair the wheel
delocate-wheel --require-archs {delocate_archs} -w {dest_dir} -v {wheel}