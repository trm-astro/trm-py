# Python Package trm-py

This is the python wrapper around trm-subs and some other trm c++ software.

This is the update to many of the repositories found [here](https://github.com/trmrsh?tab=repositories&q=&type=&language=python&sort=).
This code provides drop in replacements for the python codes.
The following four main libraries are available:

- trm_py.subs
- trm_py.roche
- trm_py.doppler
- trm_py.observing

In addition direct access to the python bindings to the cpp code is available:

- trm_py._cpp._cpp_subs
- trm_py._cpp._cpp_roche
- trm_py._cpp.cpp_doppler

This code will be installable via `pip install trm-py` from PyPi (and potentially via Conda). This code *will not* require additional C++ library installations separate from the python build system.

Until further notice this is a WIP and should not be considered usable for scientific use.

## Cibuildwheel 

> Get the main repo
> `git clone git@github.com:trm-astro/trm-py.git`
>
> Then get all the submodules
> `git submodule update --init --recursive`
>
> Install cibuildwheel 
> `pip3 install cibuildwheel`
>
> >If you need to install/start docker
> >
> > Install Docker
> > yum/brew docker (or however you like, note it must be REAL docker not podman )
> > for ubuntu/debian: `curl -sSL https://get.docker.com/ | sudo sh`
> >
> > Start Docker Deamon
> > sudo systemctl start docker
> >
> > Make sure docker works (on my vm, vagrant, I did), ymmv
> >
> > sudo groupadd docker
> > sudo usermod -aG docker $USER
> > reboot
> > 
>
> cibuildwheel


## CPP-Only build

By necessity the CPP subdirectorys are included to build the python wheel against.
This makes this the a convenient place to do a full CPP build for those looking to access trm-subs and associated libs via cpp code.
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

## Troubleshooting/Probable FAQs

PLPlot is complaining about fonts (mac):
This is an issue with plplot not bringing its dependencies with it: current fix

 > `brew install plplot`

 This brings it to the system and then it can find them. Not ideal but a good patch.

My cpp libs (conan) wont install correctly:

> 1. `cd py_subs`
> 2. `rm -r build` (potentially not required)
> 3. `conan install . --build=missing`
> 4. `cd ..`
> 5. `poetry build`

Conan is complaining about profiles

> 1. `which conan` ensure this finds conan (and where you expect it to be)
> 2. `conan profile detect` Make sure the output matches your system configuration

Conan wont find fttw3

> This seems to be a (temporary?) issue with the new conan2 repo not having the build for mac
> We can install from the old remote by adding it
>
> 1. `conan remote add old-conan https://center.conan.io`
> 2. `conan install . -r old-conan --build=missing`
>
> Now its built and cached you should have no issues

CMake 3.5 issues

> The latest CMake >4.00 deprecated anything with cmake_minimum_required(3.5)
> Install CMake 3.31
