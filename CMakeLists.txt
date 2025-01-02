# This Cmake file is used to for the pybind11 bindings of the trm librarys

# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.23)

project(trm-py)

# Find the module devlopment requirements
set(PYBIND11_FINDPYTHON ON) # Removes the need for FindPython 
find_package(pybind11 CONFIG REQUIRED)

# Target the CPP code
pybind11_add_module(trm_py src/main.cc)

install(TARGETS trm_py DESTINATION .)