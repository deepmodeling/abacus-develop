###############################################################################
# - Find cereal
# Find the native cereal headers.
#
#  CEREAL_FOUND - True if cereal is found.
#  CEREAL_INCLUDE_DIR - Where to find cereal headers.

find_path(CEREAL_INCLUDE_DIR
    cereal/cereal.hpp
    HINTS ${CEREAL_INCLUDE_DIR}
    HINTS ${Cereal_INCLUDE_DIR}
)

if(NOT CEREAL_INCLUDE_DIR)
    include(FetchContent)
    FetchContent_Declare(
        cereal
        URL https://github.com/USCiLab/cereal/archive/refs/tags/v1.3.0.tar.gz
    )
    FetchContent_Populate(cereal)
    set(CEREAL_INCLUDE_DIR ${cereal_SOURCE_DIR}/include)
endif()
# Handle the QUIET and REQUIRED arguments and
# set Cereal_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cereal DEFAULT_MSG CEREAL_INCLUDE_DIR)

# Copy the results to the output variables and target.
mark_as_advanced(CEREAL_INCLUDE_DIR)
