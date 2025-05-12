# This Cmake file is used to for the pybind11 bindings of the trm librarys

# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.27)

message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message("CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message("CMAKE_TOOlCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}")

project(py-subs)
message("cmake lists for trm_py")

#Build PLPLOT if not already built
if(NOT DEFINED CMAKE_COMMON)
    set(CMAKE_COMMON ./cmake_common)
    file(REAL_PATH "${CMAKE_COMMON}" CMAKE_COMMON)
    message("Set CMAKE_COMMON: ${CMAKE_COMMON}")
endif()
if(NOT DEFINED PLPLOT_FOUND)
    include("${CMAKE_COMMON}/plplot_install.cmake")
endif()

# The subdirectories are the libraries that are built
# The .so [or .dylib if on mac] files are copied to the lib directory within the wheel

# e.g.
# trm_py*.whl/
# ├── include/
# │   └── *.h
# ├── lib/
# │   ├── cmake/
# │   ├── libbinary.so
# │   ├── libroche.so
# │   ├── libsubs.so
# │   └── libtrm-mem.so
# ├── trm_py/
# ├── trm_py-*.dist-info/

# These must be in order to build the subdirectories
# roche -> binary -> subs
add_subdirectory(src/cpp-subs)
add_subdirectory(src/cpp-binary)
add_subdirectory(src/cpp-roche)
add_subdirectory(src/cpp-mem)



# Find the module development requirements
set(PYBIND11_FINDPYTHON ON)
message("${PYTHON_SITE_PACKAGES}")
set(pybind11_DIR "${PYTHON_SITE_PACKAGES}/pybind11/share/cmake/pybind11")
message("${pybind11_DIR}")
find_package(pybind11 CONFIG REQUIRED)
find_package(FFTW3 REQUIRED)


# Get all the targeted CPP code files
file(GLOB SUBS_SOURCES
    src/py_subs.cpp
    src/py_subs.h
)
file(GLOB ROCHE_SOURCES
    src/py_roche.cpp
    src/py_roche.h
)
file(GLOB DOPPLER_SOURCES
    src/doppler.cpp
    src/doppler.h
    src/py_doppler.cpp
    src/py_doppler.h
)

# Target the CPP code
pybind11_add_module(_cpp_subs ${SUBS_SOURCES})
target_link_libraries(_cpp_subs 
    PUBLIC
        subs
        pybind11::headers
)

pybind11_add_module(_cpp_roche ${ROCHE_SOURCES})
target_link_libraries(_cpp_roche 
    PUBLIC
        subs
        roche
        pybind11::headers
)

pybind11_add_module(_cpp_doppler ${DOPPLER_SOURCES})
target_link_libraries(_cpp_doppler 
    PUBLIC
        FFTW3::fftw3
        trm-mem
        pybind11::headers
)

target_include_directories(_cpp_subs
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_include_directories(_cpp_roche
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_include_directories(_cpp_doppler
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)


set_target_properties(_cpp_subs _cpp_roche _cpp_doppler PROPERTIES
    INSTALL_RPATH "@loader_path/../../lib"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    BUILD_RPATH_USE_ORIGIN TRUE
)


install(TARGETS _cpp_subs _cpp_roche _cpp_doppler 
    LIBRARY DESTINATION trm_py/_cpp/
    RUNTIME DESTINATION trm_py/_cpp/
    )