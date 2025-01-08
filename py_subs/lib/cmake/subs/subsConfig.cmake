

include(CMakeFindDependencyMacro)

# Import targets from the installed export set
include("${CMAKE_CURRENT_LIST_DIR}/subsTargets.cmake")

# Set any additional variables if needed, such as version or paths
set(subs_VERSION "")
