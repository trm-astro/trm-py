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

# Get the source code

> Get the main repo
> `git clone git@github.com:trm-astro/trm-py.git`
>
> `cd trm-py`
>
> Switch to your branch e.g.
>
> `git checkout <branch-name>`
>
> Then get all the submodules
>
> `git submodule update --init --recursive`

## Local Wheel Build

> Get the source as above
>
> Make a virtual environment
>
> `python3 -m venv .venv`
>
> `source .venv/bin/activate`
>
> Install scikit-build-core-conan and build
>
> `pip install scikit-build-core-conan build`
>
> `python -m build`
>
> Your wheel is now in `<trm-py>/dist/*.whl`, you will need to audit (linux) or delocate (mac) the wheel to make it portable to other systems.
>
> Note, on linux, the user running the build must have write access to `/usr/lib64` as this is the install location for the .so files it is possible to change this but requires modifying many of the CMake files looking for where LOADER_PATH is set. This will likely cause issues with the cibuildwheel so don't push this change.

## Cibuildwheel

> Get the source as above
>
> Make a virtual environment
>
> `python3 -m venv .venv`
>
> `source .venv/bin/activate` Install cibuildwheel 
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
> > `sudo groupadd docker`
> > `sudo usermod -aG docker $USER`
> > `reboot`
>
> `cibuildwheel`
>
> All your wheels are in `/wheelhouse` and in theory portable.

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

I just want the C++

> Look in WHEEL_STRUCTURES.md you will see that all the cpp is bundled in the wheel it should be possible to link against it.
> If you have issues make use of otool/ldd to check what the .dylib/.so are pointing at, all dependancies should be contained in the wheel.
> If you need more than this it should be possible to make minor alterations to the CMakefiles to enable it to install all the .dylib/.so files to the normal system locations.