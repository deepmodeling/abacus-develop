###############################################################################
# - Find ELPA
# Find the native ELPA headers and libraries.
#
#  ELPA_FOUND        - True if libelpa is found.
#  ELPA_LIBRARIES    - List of libraries when using libyaml
#  ELPA_INCLUDE_DIRS - Where to find ELPA headers.
#

find_path(ELPA_INCLUDE_DIR
    elpa/elpa.h
    HINTS ${ELPA_ROOT}
    PATH_SUFFIXES "include"
    )
find_library(ELPA_LIBRARY
    NAMES elpa
    HINTS ${ELPA_ROOT}
    PATH_SUFFIXES "lib"
    )

# Handle the QUIET and REQUIRED arguments and
# set ELPA_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ELPA DEFAULT_MSG ELPA_LIBRARY ELPA_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(ELPA_FOUND)
    set(ELPA_LIBRARIES ${ELPA_LIBRARY})
    set(ELPA_INCLUDE_DIRS ${ELPA_INCLUDE_DIR})

    if(NOT TARGET ELPA::ELPA)
        add_library(ELPA::ELPA UNKNOWN IMPORTED)
        set_target_properties(ELPA::ELPA PROPERTIES
           IMPORTED_LINK_INTERFACE_LANGUAGES "C"
           IMPORTED_LOCATION "${ELPA_LIBRARY}"
           INTERFACE_INCLUDE_DIRECTORIES "${ELPA_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(ELPA_INCLUDE_DIR ELPA_LIBRARY)
