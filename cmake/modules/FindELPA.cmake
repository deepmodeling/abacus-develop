###############################################################################
# - Find ELPA
# Find the native ELPA headers and libraries through pkg-config.
#

find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
  # Find preferred library corresponding with ABACUS configuration first
  if(ENABLE_OPENMP)
    pkg_search_module(ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa_openmp elpa)
  else()
    pkg_search_module(ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa elpa_openmp)
  endif()
else()
  message(FATAL_ERROR "Pkg-config is needed to get all information about the ELPA library")
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

if(${ELPA_VERSION} VERSION_LESS "2021.05.001")
  message(FATAL_ERROR "ELPA version >= 2021.05.001 is required.")
endif()

mark_as_advanced(ELPA_INCLUDE_DIR ELPA_LIBRARY)
