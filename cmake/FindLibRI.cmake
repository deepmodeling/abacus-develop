###############################################################################
# - Find LibRI
# Find the native LibRI files.
#
#  LIBRI_FOUND - True if LibRI is found.
#  LIBRI_DIR - Where to find LibRI files.

find_path(LIBRI_DIR
    include/RI/version.h
    HINTS ${LIBRI_DIR}
    HINTS ${LibRI_DIR}
    HINTS ${libri_DIR}
)

if(NOT LIBRI_DIR)
    include(FetchContent)
    FetchContent_Declare(
        LibRI
        URL https://github.com/abacusmodeling/LibRI/archive/refs/tags/v0.2.1.1.tar.gz
    )
    FetchContent_Populate(LibRI)
    set(LIBRI_DIR ${libri_SOURCE_DIR})
endif()
# Handle the QUIET and REQUIRED arguments and
# set LIBRI_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibRI DEFAULT_MSG LIBRI_DIR)

# Copy the results to the output variables and target.
mark_as_advanced(LIBRI_DIR)
