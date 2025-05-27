# This makes and installs the CPP executables and the shared librarys
cmake_minimum_required(VERSION 3.23)
project(trm-cpp)

include(ExternalProject)

message("CMakeLists.txt for trm-cpp")

# Determine platform-specific loader path for shared libraries
if(APPLE)
    set(LOADER_PATH "@loader_path/../../lib")
elseif(UNIX)
    # For linux we need to set the loader path to /usr/lib64
    set(LOADER_PATH "/usr/lib64")
else()
    set(LOADER_PATH "")  # Windows or unsupported platform
endif()
message("LOADER_PATH set: ${LOADER_PATH}\n")


# Define CMAKE_COMMON if not already defined
if(NOT DEFINED CMAKE_COMMON)
    set(CMAKE_COMMON ${PROJECT_SOURCE_DIR}/../cmake_common)
endif()
if(NOT DEFINED PLPLOT_FOUND)
    include("${CMAKE_COMMON}/plplot_install.cmake")
endif()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cpp-subs)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cpp-binary)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cpp-colly)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cpp-roche)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cpp-mem)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve)

# Make the lcurve executable
add_executable(lprofile ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/lprofile.cc)
add_executable(levmarq ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/levmarq.cc)
add_executable(lroche ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/lroche.cc)
#add_executable(lroches ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/lroche.cc)  # Assuming lroches is similar to lroche
add_executable(picture ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/picture.cc ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/pos_disc.cc ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/disc_eclipse.cc)
add_executable(rotprof ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/rotprof.cc)
add_executable(simplex ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/simplex.cc)
add_executable(visualise ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/src/visualise.cc)

foreach(EXEC lprofile levmarq lroche picture rotprof simplex visualise)
    #target_compile_options(${EXEC} PRIVATE ${OpenMP_CXX_FLAGS})
    target_include_directories(
        ${EXEC} PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/cpp-lcurve/include
    )

    target_link_libraries(
        ${EXEC}
            subs
            binary
            colly
            roche
            lcurve
    )

    set_target_properties(${EXEC} PROPERTIES
        LINK_FLAGS "-Wl,-rpath,${LOADER_PATH}"
        INSTALL_RPATH "${LOADER_PATH}"
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH_USE_LINK_PATH TRUE
        BUILD_RPATH_USE_ORIGIN TRUE
    )
    install(TARGETS ${EXEC} 
        DESTINATION trm_py/_cpp/
    )
    # Debugging: Print RPATH and linked libraries
    get_target_property(RPATH ${EXEC} INSTALL_RPATH)
    message("RPATH for ${EXEC}: ${RPATH}")
endforeach()

message("CMakeLists.txt for trm-cpp complete\n")