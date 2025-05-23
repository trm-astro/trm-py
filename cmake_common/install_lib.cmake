# In this file, we define the install_lib function that is used to install
# the library and its headers. The function takes the following parameter:
# - PACKAGE_NAME: The name of the package

# This is minimal as the bulk of the work is done in add_subdirectory that are called before this

# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.27)

# Function to install the library and its headers
function(install_lib PACKAGE_NAME)

    # Specify the installation directory for the library
    # For macOS this has been fine to be in lib in the wheel
    # For Linux this needs to be in lib64 in the usual system locations
    # Determine platform-specific loader path for shared libraries
    if(APPLE)
        set(LIB_PATH "lib")
            # Apply correct local linking this block assume that all the libraries are in the same directory
            set_target_properties(${PACKAGE_NAME} PROPERTIES
            INSTALL_RPATH "@loader_path"
            BUILD_WITH_INSTALL_RPATH TRUE
            INSTALL_RPATH_USE_LINK_PATH TRUE
        )
    elseif(UNIX)
        set(LIB_PATH "/usr/lib64")
    else()
        set(LIB_PATH "")  # Windows or unsupported platform
    endif()
    message("LIB_PATH set: ${LIB_PATH}\n")

    # Specify where the library binary will go
    install(TARGETS ${PACKAGE_NAME}
        EXPORT ${PACKAGE_NAME}Targets              # Export the library target for package config
        LIBRARY DESTINATION "${LIB_PATH}"              # For shared libraries
        ARCHIVE DESTINATION "${LIB_PATH}"             # For static libraries
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
        DESTINATION "${LIB_PATH}/cmake/${PACKAGE_NAME}"
    )

    # # Export the library target so other projects can import it
    install(EXPORT ${PACKAGE_NAME}Targets
        FILE ${PACKAGE_NAME}Targets.cmake
        NAMESPACE ${PACKAGE_NAME}::
        DESTINATION "${LIB_PATH}/cmake/${PACKAGE_NAME}"
    )
endfunction()

