# Wheel Structures

This file is intended to give developers a better understanding of the outcomes of the wheels built.

This is particularly important to understand with respect to the correct linking which differs from macos to linux.

In the top level CMakeLists.txt:

```CMAKE
# Determine platform-specific loader path for shared libraries
if(APPLE)
    set(LOADER_PATH "@loader_path/../../lib")
elseif(UNIX)
    # For linux we need to set the loader path to /usr/lib64
    set(LOADER_PATH "/usr/lib64")
else()
    set(LOADER_PATH "")  # Windows or unsupported platform
endif()
message("LOADER_PATH set: ${LOADER_PATH}\n")
```

Note that for mac the path is set to within the wheel on the loader path. For linux it's set to system.
`auditwheel` and `delocate` act differently meaning we can ship the `.dylib` top level in the package which is beneficial if a mac user wanted to be able to use these directly.
`auditwheel` requires these to be in a system location, we choose `/usr/lib64` then it copies them with the other system packages to the wheel.

## MacOS wheel structure

The wheel remains mostly unchanged during the delocate-wheel command in `cibw_scripts/macosx_rw.sh` the libs inspected are `_cpp_*.dylib` and the file created is `trm_py/.dylibs`.

``` bash
trm_py*.whl/
├── include/
│   └── *.h
├── lib/
│   ├── cmake/
│   ├── libbinary.dylib
│   ├── libroche.dylib
│   ├── libsubs.dylib
│   └── libtrm-mem.dylib
├── trm_py/
│   ├── _cpp/
│   │   ├── _cpp_subs.dylib
│   │   ├── _cpp_roche.dylib
│   │   ├── _cpp_doppler.dylib
│   │   ├── levmarq
│   │   ├── lprofile
│   │   ├── lroche
│   │   ├── lroches
│   │   ├── picture
│   │   ├── rotprof
│   │   ├── simplex
│   │   └── visualise
│   ├── __init__.py
│   ├── doppler/
│   │   └── <doppler submodule>
│   ├── lcurve/
│   │   └── <lcurve submodule>
│   ├── ovserving/
│   │   └── <observing submodule>
│   ├── roche/
│   │   └── <roche submodule>
│   ├── subs/
│   │   └── <subs submodule>
│   └── .dylibs/
│       └── <systemlibs>.dylib (from the wheel fix)
└── trm_py-*.dist-info/
    └── <standard python metadata>
```

## Linux Wheel Structure

We note specifically the `libbinary`, `libroche`, `libsubs`, and `libtrm-mem` that are copied from `/usr/lib64` during `auditwheel repair`

``` bash
trm_py*.whl/
├── include/
│   └── *.h
├── trm_py.libs/
│   ├── <systemlibs>.so (from the wheel fix)
│   ├── libbinary-<some_uid>.so
│   ├── libroche-<some_uid>.so
│   ├── libsubs-<some_uid>.so
│   └── libtrm-mem-<some_uid>.so
├── trm_py/
│   ├── _cpp/
│   │   ├── _cpp_subs.so
│   │   ├── _cpp_roche.so
│   │   ├── _cpp_doppler.so
│   │   ├── levmarq
│   │   ├── lprofile
│   │   ├── lroche
│   │   ├── lroches
│   │   ├── picture
│   │   ├── rotprof
│   │   ├── simplex
│   │   └── visualise
│   ├── __init__.py
│   ├── doppler/
│   │   └── <doppler submodule>
│   ├── lcurve/
│   │   └── <lcurve submodule>
│   ├── ovserving/
│   │   └── <observing submodule>
│   ├── roche/
│   │   └── <roche submodule>
│   └── subs/
│       └── <subs submodule>
└── trm_py-*.dist-info/
    └── <standard python metadata>
```

