# - Find LAPACK
# Find the native double precision ScaLAPACK headers and libraries.
#
#  LAPACK_LIBRARIES     - List of libraries when using ScaLAPACK.
#  LAPACK_FOUND         - True if ScaLAPACK is found.
#

find_library(LAPACK_LIBRARY
    NAMES openblas blas
    HINTS ${LAPACK_DIR}
    PATH_SUFFIXES "lib"
)

# Handle the QUIET and REQUIRED arguments and
# set ScaLAPACK_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK DEFAULT_MSG LAPACK_LIBRARY)

# Copy the results to the output variables and target.
if(LAPACK_FOUND)
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})

    if(NOT TARGET LAPACK::LAPACK)
        add_library(LAPACK::LAPACK UNKNOWN IMPORTED)
        set_target_properties(LAPACK::LAPACK PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${LAPACK_LIBRARY}")
    endif()
endif()

mark_as_advanced(LAPACK_LIBRARY)
