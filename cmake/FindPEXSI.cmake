###############################################################################
# - Find PEXSI
# Find PEXSI and its dependencies.
#
#  PEXSI_FOUND          - True if pexsi is found.
#  PEXSI_INCLUDE_DIR    - Where to find pexsi headers.
#  PEXSI_LIBRARY        - pexsi library.
#  ParMETIS_INCLUDE_DIR - Where to find pexsi headers.
#  ParMETIS_LIBRARY     - parmetis library.
#  METIS_LIBRARY        - metis library.
#  SuperLU_DIST_LIBRARY - superlu_dist library.

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

find_library(SuperLU_DIST_LIBRARY
    NAMES superlu_dist
    HINTS ${SuperLU_DIST_DIR}
    PATH_SUFFIXES "lib"
)

# Handle the QUIET and REQUIRED arguments and
# set PEXSI_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PEXSI DEFAULT_MSG PEXSI_LIBRARY PEXSI_INCLUDE_DIR ParMETIS_LIBRARY METIS_LIBRARY SuperLU_DIST_LIBRARY)


# Copy the results to the output variables and target.
mark_as_advanced(PEXSI_LIBRARY PEXSI_INCLUDE_DIR ParMETIS_LIBRARY SuperLU_DIST_LIBRARY)

