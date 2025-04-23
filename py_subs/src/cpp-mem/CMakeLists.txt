cmake_minimum_required(VERSION 3.23)

# Project name and version
project(cpp-mem)

include(ExternalProject)
message("CMakeLists.txt for cpp-mem")
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message("CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

set(PACKAGE_NAME mem)

# Add subdirectories for source and include files
# Add source files
set(SOURCES
    src/memsys.cc
)

# Add stubs for undefined functions
list(APPEND SOURCES src/stubs.cc)

add_library(${PACKAGE_NAME} SHARED 
    ${SOURCES}
)

target_include_directories(${PACKAGE_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)


###########################################################
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

configure_file(memConfig.cmake.in
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
    NAMESPACE mem::
    DESTINATION lib/cmake/${PACKAGE_NAME}
)