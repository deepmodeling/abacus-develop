###############################################################################
# - Find NEP_CPU
# Finds the NEP_CPU header and library.
#
# This module will search for the NEP_CPU library, looking for a hint
# from the NEP_DIR environment variable or CMake variable.
#
# This module defines the following variables:
#
#  NEP_FOUND        - True if the NEP_CPU library and headers were found.
#  NEP_INCLUDE_DIR  - The directory where nep.h is located.
#  NEP_LIBRARY      - The full path to the NEP_CPU library.
#
# It also defines the following imported target:
#
#  NEP::nep_cpu     - The NEP_CPU library target.
#
###############################################################################

find_path(NEP_INCLUDE_DIR nep.h
    HINTS ${NEP_DIR}
    PATH_SUFFIXES "include"
    )

find_library(NEP_LIBRARY
    NAMES nepcpu
    HINTS ${NEP_DIR}
    PATH_SUFFIXES "lib"
    )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NEP
    DEFAULT_MSG
    NEP_LIBRARY NEP_INCLUDE_DIR)

if(NEP_FOUND)
    if(NOT TARGET NEP::nep_cpu)
        add_library(NEP::nep_cpu UNKNOWN IMPORTED)
        set_target_properties(NEP::nep_cpu PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${NEP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${NEP_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(NEP_INCLUDE_DIR NEP_LIBRARY)