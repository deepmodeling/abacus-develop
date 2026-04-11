# - Find ScaLAPACK
# Find the native double precision ScaLAPACK headers and libraries.
#
#  ScaLAPACK_LIBRARIES     - List of libraries when using ScaLAPACK.
#  ScaLAPACK_FOUND         - True if ScaLAPACK is found.
#

# Accept common root hints from cache vars and environment.
set(_scalapack_hints
    ${SCALAPACK_DIR}
    ${SCALAPACK_ROOT}
    $ENV{SCALAPACK_DIR}
    $ENV{SCALAPACK_ROOT}
)

find_library(ScaLAPACK_LIBRARY
    NAMES
      scalapack
      scalapack-openmpi
      scalapack-mpi
      scalapack-mpich
    HINTS ${_scalapack_hints}
    PATH_SUFFIXES lib lib64
)

unset(_scalapack_hints)

# Handle the QUIET and REQUIRED arguments and
# set ScaLAPACK_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    ScaLAPACK
    REQUIRED_VARS ScaLAPACK_LIBRARY
)

# Copy the results to the output variables and target.
if(ScaLAPACK_FOUND)
    set(ScaLAPACK_LIBRARIES ${ScaLAPACK_LIBRARY})

    if(NOT TARGET ScaLAPACK::ScaLAPACK)
        add_library(ScaLAPACK::ScaLAPACK UNKNOWN IMPORTED)
        set_target_properties(ScaLAPACK::ScaLAPACK PROPERTIES
            IMPORTED_LOCATION "${ScaLAPACK_LIBRARY}")
    endif()
endif()

mark_as_advanced(ScaLAPACK_LIBRARY)
