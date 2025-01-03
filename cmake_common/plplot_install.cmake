# Configure the PLPLOT library

# Default behaviour is to download and build the PLPLOT library for the system
# PLPLOT_BUILD_TYPE:
#   Nobuild: 0 (PLPLOT_USE_PATH must be set) get lib from location
#   Local: 1 (Default) Build in the build directory
#   Path: 2 (PLPLOT_USE_PATH must be set)
#   System: 3 Build in /usr/local/src/plplot
#   Brew: 4 Use Homebrew to get the library (reccomended for mac)
# PLPLOT_USE_PATH (Required for Path or Nobuild)

include(ExternalProject)

if (NOT DEFINED PLPLOT_BUILD_TYPE)
    set(PLPLOT_BUILD_TYPE 1)
endif()

#TODO: Test this
# We need to set the name for the plplot lib
if (NOT DEFINED PLPLOT_LIB_NAME)
    set(PLPLOT_LIB_NAME libplplotcxx.dylib) # works on mac
endif()

if(NOT DEFINED PLPLOT_USE_PATH)
    if(PLPLOT_BUILD_TYPE EQUAL 0 OR PLPLOT_BUILD_TYPE EQUAL 2)
        message(FATAL_ERROR "PLPLOT_USE_PATH must be set if PLPLOT_BUILD_TYPE is ${PLPLOT_BUILD_TYPE}")
    elseif(PLPLOT_BUILD_TYPE EQUAL 1)
        set(PLPLOT_USE_PATH ${CMAKE_BINARY_DIR}/plplot_install)
    elseif(PLPLOT_BUILD_TYPE EQUAL 3)
        set(PLPLOT_USE_PATH /usr/local/src/plplot)
    endif()
endif()


if(PLPLOT_BUILD_TYPE EQUAL 0)
    message("PLPLOT_BUILD_TYPE: Nobuild")
    # Assert that PLPLOT_USE_PATH is set
    if(NOT DEFINED PLPLOT_USE_PATH)
        message(FATAL_ERROR "PLPLOT_USE_PATH must be set if PLPLOT_BUILD_TYPE is 0")
    endif()
    set(PLPLOT_LIB_PATH ${PLPLOT_USE_PATH})
elseif(PLPLOT_BUILD_TYPE EQUAL 1 OR PLPLOT_BUILD_TYPE EQUAL 2 OR PLPLOT_BUILD_TYPE EQUAL 3) # Build at defined or dfault path
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
    add_dependencies(${PACKAGE_NAME} plplot)
    set(PLPLOT_LIB_PATH ${PLPLOT_USE_PATH})

elseif(PLPLOT_BUILD_TYPE EQUAL 4) # Brew build of PLPLOT
    message("PLPLOT_BUILD_TYPE: Brew")
    # Set PLplot paths using Homebrew's default install location
    set(PLPLOT_LIB_PATH /opt/homebrew/opt/plplot)
endif()

if(NOT DEFINED PLPLOT_LIB_PATH)
    message(FATAL_ERROR "PLPLOT_LIB_PATH undefined")
else()
    set(PLPLOT_INCLUDE_DIR ${PLPLOT_LIB_PATH}/include/plplot)
    set(PLPLOT_LIB_DIR ${PLPLOT_LIB_PATH}/lib)
    set(PLPLOT_LIB ${PLPLOT_LIB_DIR}/${PLPLOT_LIB_NAME})
    set(PLPLOT_FOUND TRUE)
    message("PLPLOT_INCLUDE_DIR: ${PLPLOT_INCLUDE_DIR}")
    message("PLPLOT_LIB_DIR: ${PLPLOT_LIB_DIR}")
    message("PLPLOT_LIB: ${PLPLOT_LIB}")
endif()