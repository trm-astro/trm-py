cmake_minimum_required(VERSION 3.23)
project(py-subs)

include(ExternalProject)

message("CMakeLists.txt for py-subs")

set(PACKAGE_NAME pysubs)

# Define CMAKE_COMMON if not already defined

if(NOT DEFINED CMAKE_COMMON)
    set(CMAKE_COMMON ${PROJECT_SOURCE_DIR}/../cmake_common)
endif()
if(NOT DEFINED PLPLOT_FOUND)
    include("${CMAKE_COMMON}/plplot_install.cmake")
endif()

add_subdirectory(cpp-subs)
add_subdirectory(cpp-binary)
add_subdirectory(cpp-colly)
add_subdirectory(cpp-roche)
add_subdirectory(cpp-lcurve)

# Define sources

set(SOURCES test.cpp)

# Define the library
add_library(${PACKAGE_NAME} SHARED ${SOURCES})

# Define the include directories
target_include_directories(${PACKAGE_NAME} 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(${PACKAGE_NAME} 
    PUBLIC
        binary
        colly
        roche
        lcurve
)

# Set RPATH properties
set_target_properties(${PACKAGE_NAME} PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    BUILD_RPATH_USE_ORIGIN TRUE
)


# Make the executable
add_executable(test test.cpp)

target_link_libraries(test
    PUBLIC
        ${PACKAGE_NAME}

)

# Make the lcurve executable
add_executable(lprofile cpp-lcurve/src/lprofile.cc)
add_executable(levmarq cpp-lcurve/src/levmarq.cc)
add_executable(lroche cpp-lcurve/src/lroche.cc)
add_executable(lroches cpp-lcurve/src/lroche.cc)  # Assuming lroches is similar to lroche
add_executable(picture cpp-lcurve/src/picture.cc cpp-lcurve/src/pos_disc.cc cpp-lcurve/src/disc_eclipse.cc)
add_executable(rotprof cpp-lcurve/src/rotprof.cc)
add_executable(simplex cpp-lcurve/src/simplex.cc)
add_executable(visualise cpp-lcurve/src/visualise.cc)

foreach(EXEC lprofile levmarq lroche lroches picture rotprof simplex visualise)
    target_link_libraries(
        ${EXEC}
            subs
            binary
            colly
            roche
            lcurve
    )
endforeach()


