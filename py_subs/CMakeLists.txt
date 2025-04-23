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
add_subdirectory(src/cpp-mem)

# Find the module devlopment requirements
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


message("PY_SUBS_SOURCES: ${PY_SUBS_SOURCES}")

# Target the CPP code
pybind11_add_module(cpp_subs ${SUBS_SOURCES})
target_link_libraries(cpp_subs 
    PUBLIC
        subs
        pybind11::headers
)

pybind11_add_module(cpp_roche ${ROCHE_SOURCES})
target_link_libraries(cpp_roche 
    PUBLIC
        subs
        roche
        pybind11::headers
)

pybind11_add_module(cpp_doppler ${DOPPLER_SOURCES})
target_link_libraries(cpp_doppler 
    PUBLIC
        FFTW3::fftw3
        mem
        pybind11::headers
)

target_include_directories(cpp_subs
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_include_directories(cpp_roche
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_include_directories(cpp_doppler
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# set RPATH to PROJECT_SOURCE_DIR/install/lib and PROJECT_SOURCE_DIR/install/include
set_target_properties(cpp_subs cpp_roche cpp_doppler PROPERTIES
    INSTALL_RPATH "${PROJECT_SOURCE_DIR}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    BUILD_RPATH_USE_ORIGIN TRUE
)


install(TARGETS cpp_subs cpp_roche cpp_doppler LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/)