###############################################################################
# - Find ELPA
# Find the native ELPA headers and libraries.
#
#  ELPA_FOUND        - True if libelpa is found.
#  ELPA_LIBRARIES    - List of libraries when using libyaml
#  ELPA_INCLUDE_DIR - Where to find ELPA headers.
#

#find_path(ELPA_INCLUDE_DIR
#    elpa/elpa.h
#    HINTS ${ELPA_DIR}
#    PATH_SUFFIXES "include" "include/elpa"
#    )

find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
	pkg_search_module(Libxc REQUIRED IMPORTED_TARGET GLOBAL libxc)
else()
  message(
	  "LibXC : We need pkg-config to get all information about the libxc library")
endif()

# Handle the QUIET and REQUIRED arguments and
# set ELPA_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libxc DEFAULT_MSG Libxc_LINK_LIBRARIES Libxc_INCLUDE_DIRS)

# Copy the results to the output variables and target.
if(Libxc_FOUND AND NOT TARGET Libxc::xc)
	set(Libxc_LIBRARY ${Libxc_LINK_LIBRARIES})
	set(Libxc_LIBRARIES ${Libxc_LIBRARY})
	set(Libxc_INCLUDE_DIR ${Libxc_INCLUDE_DIRS})
	add_library(Libxc::xc UNKNOWN IMPORTED)
	set_target_properties(Libxc::xc PROPERTIES
		IMPORTED_LOCATION "${Libxc_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES "${Libxc_INCLUDE_DIR}")
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${Libxc_INCLUDE_DIR})

mark_as_advanced(Libxc_INCLUDE_DIR Libxc_LIBRARY)
