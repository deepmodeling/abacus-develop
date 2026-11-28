###############################################################################
# - Find ELPA
# Find the native ELPA headers and libraries through pkg-config.
#

find_package(PkgConfig REQUIRED)

# Find preferred library corresponding with ABACUS configuration first
if(ENABLE_OPENMP)
  pkg_search_module(ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa_openmp elpa)
else()
  pkg_search_module(ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa)
endif()

if(ELPA_VERSION VERSION_LESS "2021.05.001")
  message(FATAL_ERROR "ELPA version >= 2021.05.001 is required.")
endif()

# Handle the QUIET and REQUIRED arguments and
# set ELPA_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ELPA DEFAULT_MSG ELPA_LINK_LIBRARIES
                                  ELPA_INCLUDE_DIRS)

if(ELPA_FOUND AND NOT TARGET ELPA::ELPA)
  add_library(ELPA::ELPA INTERFACE IMPORTED)
  # TODO: Clean up INTERFACE_COMPILE_OPTIONS
  set_target_properties(
    ELPA::ELPA
    PROPERTIES INTERFACE_LINK_LIBRARIES PkgConfig::ELPA
               INTERFACE_INCLUDE_DIRECTORIES "${ELPA_INCLUDE_DIRS}"
               INTERFACE_COMPILE_OPTIONS "${ELPA_CFLAGS_OTHER}")
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ELPA_INCLUDE_DIRS})

mark_as_advanced(ELPA_INCLUDE_DIR ELPA_LIBRARY)
