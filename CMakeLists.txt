# This Cmake file is used to for the pybind11 bindings of the trm librarys

# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.27)

project(trm-py)

add_subdirectory(src/cpp-subs)

# Find the module devlopment requirements
set(PYBIND11_FINDPYTHON ON)
set(pybind11_DIR "${PROJECT_SOURCE_DIR}/.venv/lib/python3.13/site-packages/pybind11/share/cmake/pybind11")
message("${pybind11_DIR}")
find_package(pybind11 CONFIG REQUIRED)
# Target the CPP code
pybind11_add_module(trm_py "${PROJECT_SOURCE_DIR}/src/main.cpp")
target_link_libraries(
    trm_py 
    PRIVATE 
        subs
        pybind11::headers
)

install(TARGETS trm_py DESTINATION .)