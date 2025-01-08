#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "subs::subs" for configuration "Release"
set_property(TARGET subs::subs APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(subs::subs PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsubs.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libsubs.dylib"
  )

list(APPEND _cmake_import_check_targets subs::subs )
list(APPEND _cmake_import_check_files_for_subs::subs "${_IMPORT_PREFIX}/lib/libsubs.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
