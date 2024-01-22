###############################################################################
# - Find ELPA
# Find the native ELPA headers and libraries.
#
#  ELPA_FOUND        - True if libelpa is found.
#  ELPA_LIBRARIES    - List of libraries when using libyaml
#  ELPA_INCLUDE_DIR - Where to find ELPA headers.
#

find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
  if(DEFINED ELPA_DIR)
    string(APPEND CMAKE_PREFIX_PATH ";${ELPA_DIR}")
  endif()
  if(USE_OPENMP)
    pkg_search_module(ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa_openmp)
  else()
    pkg_search_module(ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa)
  endif()
else()
  find_path(ELPA_INCLUDE_DIRS
    elpa/elpa.h
    HINTS ${ELPA_DIR}
    PATH_SUFFIXES "include" "include/elpa"
    )
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
  #message(
  #  "ELPA : We need pkg-config to get all information about the elpa library")
endif()

# Handle the QUIET and REQUIRED arguments and
# set ELPA_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ELPA DEFAULT_MSG ELPA_LINK_LIBRARIES ELPA_INCLUDE_DIRS)

# Copy the results to the output variables and target.
if(ELPA_FOUND)
    set(ELPA_LIBRARY ${ELPA_LINK_LIBRARIES})
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
