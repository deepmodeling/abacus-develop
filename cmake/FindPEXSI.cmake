###############################################################################
# - Find cereal
# Find the native cereal headers.
#
#  PEXSI_FOUND - True if cereal is found.
#  PEXSI_INCLUDE_DIR - Where to find cereal headers.

find_path(PEXSI_INCLUDE_DIR
    NAMES c_pexsi_interface.h
    HINTS ${PEXSI_DIR}
    PATH_SUFFIXES "include"
)

find_library(PEXSI_LIBRARY
    NAMES pexsi
    HINTS ${PEXSI_DIR}
    PATH_SUFFIXES "lib"
)

find_path(ParMETIS_INCLUDE_DIR
    NAMES metis.h parmetis.h
    HINTS ${ParMETIS_DIR}
    PATH_SUFFIXES "include"
)

find_library(METIS_LIBRARY
    NAMES metis
    HINTS ${ParMETIS_DIR}
    PATH_SUFFIXES "lib"
)

find_library(ParMETIS_LIBRARY
    NAMES parmetis
    HINTS ${ParMETIS_DIR}
    PATH_SUFFIXES "lib"
)

find_library(SuperLU_LIBRARY
    NAMES libsuperlu_dist.a
    HINTS ${SuperLU_DIR}
    PATH_SUFFIXES "lib"
)

# Handle the QUIET and REQUIRED arguments and
# set Cereal_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PEXSI DEFAULT_MSG PEXSI_LIBRARY PEXSI_INCLUDE_DIR ParMETIS_LIBRARY METIS_LIBRARY SuperLU_LIBRARY)


# Copy the results to the output variables and target.
mark_as_advanced(PEXSI_LIBRARY PEXSI_INCLUDE_DIR ParMETIS_LIBRARY SuperLU_LIBRARY)

