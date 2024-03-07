###############################################################################
# - Find ELPA
# Find the native ELPA headers and libraries.
#
#  ELPA_FOUND        - True if libelpa is found.
#  ELPA_LIBRARIES    - List of libraries when using libyaml
#  ELPA_INCLUDE_DIR - Where to find ELPA headers.
#

find_package(PkgConfig)

# Compatible layer towards old manual routines
if(DEFINED ELPA_INCLUDE_DIR)
  set(ELPA_INCLUDE_DIRS ${ELPA_INCLUDE_DIR})
endif()
if(DEFINED ELPA_LIBRARIES)
  set(ELPA_LINK_LIBRARIES ${ELPA_LIBRARIES})
endif()

find_path(ELPA_INCLUDE_DIRS
    elpa/elpa.h
    HINTS ${ELPA_DIR}
    PATH_SUFFIXES "include" "include/elpa"
    )
# Fix #3589
# First if judges if ELPA dir specified
if(ELPA_INCLUDE_DIRS MATCHES "^/usr/include/elpa/.*")
  # Second if judges if global visible ELPA header found
  if(DEFINED ELPA_DIR OR CMAKE_PREFIX_PATH MATCHES ".*elpa.*")
    unset(ELPA_INCLUDE_DIRS)
  endif()
endif()
if(USE_OPENMP)
    find_library(ELPA_LINK_LIBRARIES
    NAMES elpa_openmp elpa
    HINTS ${ELPA_DIR}
    PATH_SUFFIXES "lib"
    )
else()
    find_library(ELPA_LINK_LIBRARIES
    NAMES elpa
    HINTS ${ELPA_DIR}
    PATH_SUFFIXES "lib"
    )
endif()

# Incompatible with ELPA earlier than 2021.11.001
# Before ELPA 2021.11.001, its pkg-config file 
# is named like "elpa-2021.05.002.pc".
if(NOT ELPA_INCLUDE_DIRS AND PKG_CONFIG_FOUND)
  if(DEFINED ELPA_DIR)
    string(APPEND CMAKE_PREFIX_PATH ";${ELPA_DIR}")
  endif()
  if(USE_OPENMP)
    pkg_search_module(ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa_openmp)
  else()
    pkg_search_module(ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa)
  endif()
elseif(NOT PKG_CONFIG_FOUND)
  message(
    "ELPA : We need pkg-config to get all information about the elpa library")
endif()

# Handle the QUIET and REQUIRED arguments and
# set ELPA_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ELPA DEFAULT_MSG ELPA_LINK_LIBRARIES ELPA_INCLUDE_DIRS)

# Copy the results to the output variables and target.
if(ELPA_FOUND)
    list(GET ELPA_LINK_LIBRARIES 0 ELPA_LIBRARY)
    set(ELPA_INCLUDE_DIR ${ELPA_INCLUDE_DIRS})

    if(NOT TARGET ELPA::ELPA)
        add_library(ELPA::ELPA UNKNOWN IMPORTED)
        set_target_properties(ELPA::ELPA PROPERTIES
           IMPORTED_LINK_INTERFACE_LANGUAGES "C"
           IMPORTED_LOCATION "${ELPA_LIBRARY}"
           INTERFACE_INCLUDE_DIRECTORIES "${ELPA_INCLUDE_DIR}")
    endif()
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ELPA_INCLUDE_DIR})

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
#include <elpa/elpa_version.h>
#if ELPA_API_VERSION < 20210430
#error ELPA version is too old.
#endif
int main(){}
"
ELPA_VERSION_SATISFIES
)
if(NOT ELPA_VERSION_SATISFIES)
    message(FATAL_ERROR "ELPA version is too old. We support version 2021 or higher.")
endif()

mark_as_advanced(ELPA_INCLUDE_DIR ELPA_LIBRARY)
