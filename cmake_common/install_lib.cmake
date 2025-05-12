# In this file, we define the install_lib function that is used to install
# the library and its headers. The function takes the following parameter:
# - PACKAGE_NAME: The name of the package

# This is minimal as the bulk of the work is done in add_subdirectory that are called before this

# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.27)

# Function to install the library and its headers
function(install_lib PACKAGE_NAME)

    # Apply correct local linking
    set_target_properties(${PACKAGE_NAME} PROPERTIES
        INSTALL_RPATH "@loader_path"
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )

    # Specify where the library binary will go
    install(TARGETS ${PACKAGE_NAME}
        EXPORT ${PACKAGE_NAME}Targets              # Export the library target for package config
        LIBRARY DESTINATION lib              # For shared libraries
        ARCHIVE DESTINATION lib              # For static libraries
        RUNTIME DESTINATION bin              # For executables (if any)
        INCLUDES DESTINATION include         # Where to put public headers
    )

    # Header files
    # Use exists to check if the include directory exists
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
        install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
            DESTINATION include
        )
    endif()

    # Create and install CMake package configuration files
    include(CMakePackageConfigHelpers)

    # # Configure the version file
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}ConfigVersion.cmake"
        VERSION 0.1
        COMPATIBILITY AnyNewerVersion
    )

    configure_file("${PACKAGE_NAME}Config.cmake.in"
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
        NAMESPACE ${PACKAGE_NAME}::
        DESTINATION lib/cmake/${PACKAGE_NAME}
    )
endfunction()

