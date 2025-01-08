# Python Package trm-py

This is the python wrapper around trm-subs and some other trm c++ software.

This is the update to many of the repositories found [here](https://github.com/trmrsh?tab=repositories&q=&type=&language=python&sort=). The code will create a drop in replacement via the form `from trm_py import trm_<name>` or `from trm_py.trm_<name> import <name>`, `import trm_py` will also be supported with tool access via `trm_py.<name>.<tool>`.

This code will be installable via `pip install trm-py` from PyPi (and potentially via Conda). This code will not require additional C++ libaray installations as it will come formost as a `bdist` for most common platforms and a `sdist` to support build and installation on uncommon platforms.

Until further notice this is a WIP and should not be considered usable for scientific use.

## CPP Subdirectorys

By necessity the CPP subdirectorys are included to build the python wheel against.
This makes this the aconvenient place to do a full CPP build for those looking to access trm-subs and associated libs via cpp code.
To build (and install) all the TRM cpp libs directly:

TODO: Test this on multiplatform builds:

0. Decide on your PLPLOT install method (see cmake_common/plplot_install.cmake)
    1. If PLPLOT is installed then `0` and set path (recommended for linux)
    2. Build from source during the install select `1`,`2`, or `3` depending on desired build location and set path if required
    3. Use homebrew install `4` (recommended for mac) (`brew install plplot`)
1. Install/activate the python envronment using poetry (required for conan)
2. `cd src`, Move to the source directory
3. `conan install . --build=missing`, Install PCRE2, llvm-openmp, and SOFA from conan
4. `conan profile detect --force`, Set the conan build profile
5. `cmake --preset conan-release -DPLPLOT_BUILD_TYPE='0,1,2,3,or 4' -DPLPLOT_USE_PATH='path/to/plplot/install'`, use the PLPLOT build choice and use the PLPLOT path if needed
6. `cmake --build --preset conan-release`

OPTIONAL: Install to system, else programs are in `src/build/release`

7. `<sudo> cmake --install build/Release`
8. (Optional): use `otool -l <exe or lib> | grep RPATH -A2` to ensure that the programs and libraries have linked correctly
