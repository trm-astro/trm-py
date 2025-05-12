# This makes and installs the CPP executables and the shared librarys
cmake_minimum_required(VERSION 3.23)
project(trm-cpp)

include(ExternalProject)

# Define CMAKE_COMMON if not already defined
if(NOT DEFINED CMAKE_COMMON)
    set(CMAKE_COMMON ${PROJECT_SOURCE_DIR}/../cmake_common)
endif()
if(NOT DEFINED PLPLOT_FOUND)
    include("${CMAKE_COMMON}/plplot_install.cmake")
endif()
# This is a common library install script for subs, binary, colly, roche and lcurve
# The include brings in a function to install the libraries
include("${CMAKE_COMMON}/install_lib.cmake") 

add_subdirectory(cpp-subs)
add_subdirectory(cpp-binary)
add_subdirectory(cpp-colly)
add_subdirectory(cpp-roche)
add_subdirectory(cpp-lcurve)


# Use a common cmake script to intall the above to the system

foreach(DEP cpp-subs cpp-binary cpp-colly cpp-roche cpp-lcurve)
    
endforeach()


# ############################################################
# Define sources
set(SOURCES test/test.cxx)

# Define the library
add_library(${PACKAGE_NAME} SHARED ${SOURCES})

# Define the include directories
target_include_directories(${PACKAGE_NAME} 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(${PACKAGE_NAME} 
    PUBLIC
        roche
        binary
        colly
        lcurve
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
    set_target_properties(${EXEC} PROPERTIES
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH_USE_LINK_PATH TRUE
        BUILD_RPATH_USE_ORIGIN TRUE
    )
    install(TARGETS ${EXEC} DESTINATION bin)
endforeach()



