###############################################################################
# - Find cereal
# Find the native cereal headers.
#
#  CEREAL_FOUND - True if cereal is found.
#  CEREAL_INCLUDE_DIRS - Where to find cereal headers.
#

find_path(Cereal_INCLUDE_DIR
    cereal/cereal.hpp
    HINTS ${CEREAL_ROOT}
    PATH_SUFFIXES "include"
    )

# Handle the QUIET and REQUIRED arguments and
# set Cereal_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cereal DEFAULT_MSG Cereal_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(Cereal_FOUND)
    set(Cereal_INCLUDE_DIRS ${Cereal_INCLUDE_DIR})
endif()

mark_as_advanced(Cereal_INCLUDE_DIR)
