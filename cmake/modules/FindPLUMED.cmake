###############################################################################
# - Find PLUMED
# Find the native PLUMED headers and libraries.
#
# PLUMED is built with Autotools and installs a pkg-config file (plumed.pc), so
# pkg-config is the preferred discovery mechanism; a manual PLUMED_ROOT hint is
# kept as a fallback for non-standard installations.
#
# This module provides the imported target
#   PLUMED::plumed
# which carries the include directories and the wrapper library (libplumed) of
# a PLUMED installation.
#

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PLUMED_PKG QUIET IMPORTED_TARGET GLOBAL plumed)
endif()

if(PLUMED_PKG_FOUND)
  set(PLUMED_INCLUDE_DIRS ${PLUMED_PKG_INCLUDE_DIRS})
  set(PLUMED_LIBRARIES ${PLUMED_PKG_LINK_LIBRARIES})
  set(PLUMED_VERSION ${PLUMED_PKG_VERSION})
else()
  # Fallback: manual discovery, e.g. -DPLUMED_ROOT=/path/to/prefix
  find_path(PLUMED_INCLUDE_DIRS
    plumed/wrapper/Plumed.h
    HINTS ${PLUMED_ROOT}
    PATH_SUFFIXES "include"
  )
  find_library(PLUMED_LIBRARIES
    NAMES plumed
    HINTS ${PLUMED_ROOT}
    PATH_SUFFIXES "lib"
  )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PLUMED
  REQUIRED_VARS PLUMED_LIBRARIES PLUMED_INCLUDE_DIRS
  VERSION_VAR PLUMED_VERSION)

if(PLUMED_FOUND AND NOT TARGET PLUMED::plumed)
  add_library(PLUMED::plumed UNKNOWN IMPORTED)
  set_target_properties(PLUMED::plumed PROPERTIES
    IMPORTED_LOCATION "${PLUMED_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${PLUMED_INCLUDE_DIRS}"
  )
endif()

mark_as_advanced(PLUMED_INCLUDE_DIRS PLUMED_LIBRARIES)
