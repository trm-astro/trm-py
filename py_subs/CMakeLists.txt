# This Cmake file is used to for the pybind11 bindings of the trm librarys

# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.27)

project(py-subs)
message("cmake lists for py_subs")

# These must be in order to build the subdirectories
# roche -> binary -> subs
add_subdirectory(src/cpp-subs)
add_subdirectory(src/cpp-binary)
add_subdirectory(src/cpp-roche)

# Find the module devlopment requirements
set(PYBIND11_FINDPYTHON ON)
message("${PYTHON_SITE_PACKAGES}")
set(pybind11_DIR "${PYTHON_SITE_PACKAGES}/pybind11/share/cmake/pybind11")
message("${pybind11_DIR}")
find_package(pybind11 CONFIG REQUIRED)


# Get all the targeted CPP code files
file(GLOB PY_SUBS_SOURCES
    src/*.cpp
    src/*.h
)

message("PY_SUBS_SOURCES: ${PY_SUBS_SOURCES}")

# Target the CPP code
pybind11_add_module(py_subs ${PY_SUBS_SOURCES})


target_link_libraries(py_subs 
    PUBLIC
        subs
        roche
        pybind11::headers
)

target_include_directories(py_subs
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# set RPATH to PROJECT_SOURCE_DIR/install/lib and PROJECT_SOURCE_DIR/install/include
set_target_properties(py_subs PROPERTIES
    INSTALL_RPATH "${PROJECT_SOURCE_DIR}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    BUILD_RPATH_USE_ORIGIN TRUE
)

install(TARGETS py_subs LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/)