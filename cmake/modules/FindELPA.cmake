###############################################################################
# - Find ELPA
# Find the native ELPA headers and libraries through pkg-config.
#

# ========================================================================
# Deprecated (TODO: Remove this part)
# ========================================================================

# Migrate caches created when FindPkgConfig used the public ELPA prefix.
if(DEFINED CACHE{ELPA_LIBRARIES})
  get_property(_elpa_libraries_cache_type CACHE ELPA_LIBRARIES PROPERTY TYPE)
  if(_elpa_libraries_cache_type STREQUAL "INTERNAL")
    unset(ELPA_LIBRARIES CACHE)
    unset(ELPA_INCLUDE_DIRS CACHE)
    unset(ELPA_LINK_LIBRARIES CACHE)
  endif()
endif()
unset(_elpa_libraries_cache_type)

# Compatible layer towards old manual routines
if(DEFINED ELPA_DIR)
  message(WARNING "ELPA_DIR is deprecated and will be removed in the future release.")
endif()
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
if(ENABLE_OPENMP)
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

# ========================================================================

# TODO: Make pkg-config discovery unconditional after removing the deprecated
# manual discovery path.
if(NOT ELPA_INCLUDE_DIRS)
  find_package(PkgConfig)
  if(NOT PKG_CONFIG_FOUND)
    message(FATAL_ERROR "Pkg-config is needed to get all information about the ELPA library")
  endif()
  # Find preferred library corresponding with ABACUS configuration first
  if(ENABLE_OPENMP)
    pkg_search_module(ELPA_PKG REQUIRED IMPORTED_TARGET GLOBAL elpa_openmp elpa)
  else()
    pkg_search_module(ELPA_PKG REQUIRED IMPORTED_TARGET GLOBAL elpa)
  endif()
  if(${ELPA_PKG_VERSION} VERSION_LESS "2021.05.001")
    message(FATAL_ERROR "ELPA version >= 2021.05.001 is required.")
  endif()
  set(ELPA_INCLUDE_DIRS ${ELPA_PKG_INCLUDE_DIRS})
  set(ELPA_LINK_LIBRARIES ${ELPA_PKG_LINK_LIBRARIES})
  set(ELPA_LIBRARIES ${ELPA_PKG_LIBRARIES})
  set(ELPA_VERSION ${ELPA_PKG_VERSION})
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
    # TODO: Remove the manual target fallback with the deprecated ELPA inputs.
    if(TARGET PkgConfig::ELPA_PKG)
      add_library(ELPA::ELPA ALIAS PkgConfig::ELPA_PKG)
    else()
      add_library(ELPA::ELPA UNKNOWN IMPORTED)
      set_target_properties(ELPA::ELPA PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${ELPA_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${ELPA_INCLUDE_DIR}")
    endif()
  endif()
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ELPA_INCLUDE_DIR})

# Compability workaround for ELPA_DIR
# TODO: Remove this check with the deprecated ELPA_DIR path.
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
