cmake_minimum_required(VERSION 3.23)
project(cpp-subs)

include(ExternalProject)

message("CMakeLists.txt for cpp-subs")
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message("CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

set(PACKAGE_NAME subs)

############################################################
# Define sources
############################################################

set(SOURCES
    src/gauss.cc
    src/ran.cc
    src/ephem.cc
    src/gammln.cc
    src/gammp.cc
    src/gammq.cc
    src/gser.cc
    src/hitem.cc
    src/header.cc
    src/date.cc
    src/time.cc
    src/position.cc
    src/star.cc
    src/binary.cc
    src/telescope.cc
    src/svdfit.cc
    src/filnam.cc
    src/tchi.cc
    src/input.cc
    src/string_funcs.cc
    src/centroid.cc
    src/bsstep.cc
    src/mmid.cc
    src/pzextr.cc
    src/gcf.cc
    src/lud.cc
    src/gaussj.cc
    src/jacob.cc
    src/sleep.cc
    src/sigma_reject.cc
    src/what_colour.cc
    src/planck.cc
    src/fft.cc
    src/llsqr.cc
    src/fasper.cc
    src/complex.cc
    src/plot.cc
    src/formula.cc
    src/fraction.cc
    src/units.cc
    src/format.cc
    src/poly.cc
    src/rebin.cc
    src/amoeba.cc
    src/genetic.cc
    src/rtsafe.cc
    src/brent.cc
    src/dbrent.cc
    src/mnbrak.cc
    src/powell.cc
    src/safunc.cc
    src/poisson.cc
    src/extinct.cc
    src/byte_swap.cc
    src/endian.cc
    src/boxcar.cc
    src/numdiff.cc
    src/factln.cc
    src/runge_kutta.cc
    src/voigt.cc
    src/stoerm.cc
)

############################################################
# Create a library
############################################################
add_library(${PACKAGE_NAME} SHARED 
    ${SOURCES}
)

target_include_directories(${PACKAGE_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

############################################################
# Check for and install dependencies
############################################################
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/Release/generators/conan_toolchain.cmake")
    include("${CMAKE_TOOLCHAIN_FILE}")
endif()

# Configure the PLPLOT library
if(NOT DEFINED CMAKE_COMMON)
    set(CMAKE_COMMON ${PROJECT_SOURCE_DIR}/../../cmake_common)
endif()
if(NOT DEFINED PLPLOT_FOUND)
    include("${CMAKE_COMMON}/plplot_install.cmake")
endif()

find_package(SOFA REQUIRED CONFIG)
find_package(PCRE2 REQUIRED CONFIG)

# if (NOT SOFA_FOUND)
#     message(FATAL_ERROR "SOFA library not found.")
# endif()
# if (NOT PCRE2_FOUND)
#     message(FATAL_ERROR "PCRE2 library not found.")
# endif()
# if (NOT PLPLOT_FOUND)
#     message(FATAL_ERROR "PLPLOT library not found.")
# endif()

target_include_directories(${PACKAGE_NAME}
    PUBLIC 
        $<BUILD_INTERFACE:${PLPLOT_INCLUDE_DIR}>
)

target_link_libraries(${PACKAGE_NAME} 
    PUBLIC 
        sofa::sofa 
        pcre2::pcre2
        ${PLPLOT_LIB_DIR}/libplplotcxx.dylib
)

# Set RPATH properties
set_target_properties(${PACKAGE_NAME} PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    BUILD_RPATH_USE_ORIGIN TRUE
)


############################################################
# Create an executable
############################################################

# No executable in this code, this is a library


############################################################
# Install
############################################################
set_target_properties(${EXEC} PROPERTIES
INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
BUILD_WITH_INSTALL_RPATH TRUE
INSTALL_RPATH_USE_LINK_PATH TRUE
BUILD_RPATH_USE_ORIGIN TRUE
)
# Library
# Note: may not work on windows
# Specify where the library binary will go
install(TARGETS ${PACKAGE_NAME}
    EXPORT ${PACKAGE_NAME}Targets              # Export the library target for package config
    LIBRARY DESTINATION lib              # For shared libraries
    ARCHIVE DESTINATION lib              # For static libraries
    RUNTIME DESTINATION bin              # For executables (if any)
    INCLUDES DESTINATION include         # Where to put public headers
)

# Header files
install(DIRECTORY include/ 
    DESTINATION include
)

# # Create and install CMake package configuration files
include(CMakePackageConfigHelpers)

# # Configure the version file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}ConfigVersion.cmake"
    VERSION 0.1
    COMPATIBILITY AnyNewerVersion
)

configure_file(subsConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake"
    @ONLY
)

# # Install the config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}ConfigVersion.cmake"
    DESTINATION lib/cmake/${PACKAGE_NAME}
)

# # Export the library target so other projects can import it
install(EXPORT ${PACKAGE_NAME}Targets
    FILE ${PACKAGE_NAME}Targets.cmake
    NAMESPACE subs::
    DESTINATION lib/cmake/${PACKAGE_NAME}
)