# Configure the PLPLOT library

# Only 0 and 4 currently work, 1,2,3 in theory should work but are not currently working

# Default behaviour is to download and build the PLPLOT library for the system
# PLPLOT_BUILD_TYPE:
#   Nobuild: 0 (PLPLOT_USE_PATH must be set) get lib from known location
#   Local: 1 Build in the build directory
#   Submodule: 2 Use the local git submodule to build in the build directory
#   Path: 3 (PLPLOT_USE_PATH must be set)
#   PkgMgr: 4 Use package manager to install PLPLOT (Homebrew, apt, vcpkg)
# PLPLOT_USE_PATH (Required for Path or Nobuild)

include(ExternalProject)

if (NOT DEFINED PLPLOT_BUILD_TYPE)
    set(PLPLOT_BUILD_TYPE 4)
endif()

#TODO: Test this
# We need to set the name for the plplot lib
if (NOT DEFINED PLPLOT_LIB_NAME)
    if (WIN32)
        set(PLPLOT_LIB_NAME libplplotcxx.dll) 
    elseif (UNIX)
        set(PLPLOT_LIB_NAME libplplotcxx.so) 
    elseif (APPLE)
        set(PLPLOT_LIB_NAME libplplotcxx.dylib) # works on mac
    endif()
endif()

if(NOT DEFINED PLPLOT_USE_PATH)
    if(PLPLOT_BUILD_TYPE EQUAL 0 OR PLPLOT_BUILD_TYPE EQUAL 3)
        message(FATAL_ERROR "PLPLOT_USE_PATH must be set if PLPLOT_BUILD_TYPE is ${PLPLOT_BUILD_TYPE}")
    elseif(PLPLOT_BUILD_TYPE EQUAL 1 OR PLPLOT_BUILD_TYPE EQUAL 2)
        set(PLPLOT_USE_PATH ${CMAKE_BINARY_DIR}/plplot_install)
    endif()
endif()


if(PLPLOT_BUILD_TYPE EQUAL 0)
    message("PLPLOT_BUILD_TYPE: Nobuild")
    # Assert that PLPLOT_USE_PATH is set
    if(NOT DEFINED PLPLOT_USE_PATH)
        message(FATAL_ERROR "PLPLOT_USE_PATH must be set if PLPLOT_BUILD_TYPE is 0")
    endif()
    set(PLPLOT_LIB_PATH ${PLPLOT_USE_PATH})
elseif(PLPLOT_BUILD_TYPE EQUAL 2) # Build from submodule
    message("PLPLOT_BUILD_TYPE: Submodule")
    # Assert that PLPLOT_USE_PATH is set
    if(NOT DEFINED PLPLOT_USE_PATH)
        message(FATAL_ERROR "PLPLOT_USE_PATH must be set if PLPLOT_BUILD_TYPE is 2")
    endif()
    message("Building Plplot from submodule at ${CMAKE_CURRENT_SOURCE_DIR}/src/plplot")
    ExternalProject_Add(
        plplot
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/plplot
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PLPLOT_USE_PATH}
                -DCMAKE_BUILD_TYPE=Release
                -DPLPLOT_BUILD_SHARED=ON       # Optionally build shared libraries
                -DENABLE_python=OFF            # Disable Python bindings if not needed
                -DENABLE_tcl=OFF               # Disable Tcl bindings if not needed
        INSTALL_DIR ${PLPLOT_USE_PATH}
    )
    set(PLPLOT_LIB_PATH ${PLPLOT_USE_PATH})
elseif(PLPLOT_BUILD_TYPE EQUAL 1 OR PLPLOT_BUILD_TYPE EQUAL 3) # Build at defined or dfault path
    # Assert that PLPLOT_USE_PATH is set
    if(NOT DEFINED PLPLOT_USE_PATH)
        message(FATAL_ERROR "PLPLOT_USE_PATH undefined")
    else()
        message("Building Plplot at ${PLPLOT_USE_PATH}")
    endif()
    ExternalProject_Add(
        plplot
        URL https://sourceforge.net/projects/plplot/files/plplot/5.15.0%20Source/plplot-5.15.0.tar.gz
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PLPLOT_USE_PATH}
                -DCMAKE_BUILD_TYPE=Release
                -DPLPLOT_BUILD_SHARED=ON       # Optionally build shared libraries
                -DENABLE_python=OFF            # Disable Python bindings if not needed
                -DENABLE_tcl=OFF               # Disable Tcl bindings if not needed
        INSTALL_DIR ${PLPLOT_USE_PATH}
    )
    set(PLPLOT_LIB_PATH ${PLPLOT_USE_PATH})

elseif(PLPLOT_BUILD_TYPE EQUAL 4) # PkgManager build of PLPLOT
    # Get current platform
    if(APPLE)
        message("PLPLOT_BUILD_TYPE: Brew")
        # Set PLplot paths using Homebrew's default install location
        set(PLPLOT_LIB_PATH /opt/homebrew/opt/plplot)
    elseif(UNIX)
        message("PLPLOT_BUILD_TYPE: apt")
        message("Looking for: ${PLPLOT_LIB_NAME}")
        # Manually search for the plplot library
        # it's going to be in usr/lib
        # then we need the achitecture specific path
        if (EXISTS "/usr/lib/x86_64-linux-gnu/${PLPLOT_LIB_NAME}")
            set(PLPLOT_LIB_PATH /usr/lib/x86_64-linux-gnu/)
        elseif (EXISTS "/usr/lib/aarch64-linux-gnu/${PLPLOT_LIB_NAME}")
            set(PLPLOT_LIB_PATH /usr/lib/aarch64-linux-gnu/)
        elseif (EXISTS "/usr/lib/${PLPLOT_LIB_NAME}")
            set(PLPLOT_LIB_PATH /usr/lib/)
        elseif (EXISTS "/usr/local/lib/${PLPLOT_LIB_NAME}")
            set(PLPLOT_LIB_PATH /usr/local/lib/)
        else()
            message(FATAL_ERROR "PLPLOT_LIB_PATH not found")
        endif()


        #cmake_path(REMOVE_FILENAME "${PLPLOT_FULL_PATH}" OUTPUT_VARIABLE PLPLOT_LIB_PATH)
        message("cmake resolved parent as: PLPLOT_LIB_PATH: ${PLPLOT_LIB_PATH}")
    elseif(WIN32)
        message("PLPLOT_BUILD_TYPE: vcpkg")
        # Set PLplot paths using vcpkg's default install location
        find_package(plplot CONFIG REQUIRED)
    else()
        message(FATAL_ERROR "PLPLOT_BUILD_TYPE: Unsupported platform")
    endif()
endif()

if(NOT DEFINED PLPLOT_LIB_PATH)
    message(FATAL_ERROR "PLPLOT_LIB_PATH undefined")
else()
    message("pre-resolve PLPLOT_LIB_PATH: ${PLPLOT_LIB_PATH}")
    file(REAL_PATH "${PLPLOT_LIB_PATH}" PLPLOT_LIB_PATH)
    set(PLPLOT_INCLUDE_DIR ${PLPLOT_LIB_PATH}/include/plplot)
    set(PLPLOT_LIB_DIR ${PLPLOT_LIB_PATH}/lib)
    file(REAL_PATH "${PLPLOT_LIB_DIR}/${PLPLOT_LIB_NAME}" PLPLOT_LIB)
    set(PLPLOT_FOUND TRUE)
    message("PLPLOT_INCLUDE_DIR: ${PLPLOT_INCLUDE_DIR}")
    message("PLPLOT_LIB_DIR: ${PLPLOT_LIB_DIR}")
    message("PLPLOT_LIB: ${PLPLOT_LIB}")
endif()