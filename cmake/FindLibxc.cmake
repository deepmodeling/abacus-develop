include(FindPackageHandleStandardArgs)

if(DEFINED Libxc_DIR)
  string(APPEND CMAKE_PREFIX_PATH ";${Libxc_DIR}")
endif()
# Using pkg-config interface as default, to
# avoid linking to wrong global visible Libxc instead of
# specified one.
# NO REQUIRED here, otherwhile it would throw error
# with no LibXC found.
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
  pkg_search_module(Libxc IMPORTED_TARGET GLOBAL libxc)
  find_package_handle_standard_args(Libxc DEFAULT_MSG Libxc_LINK_LIBRARIES Libxc_INCLUDE_DIRS)
endif()
if(NOT Libxc_FOUND)
  find_package(Libxc REQUIRED HINTS
    ${Libxc_DIR}/share/cmake/Libxc
    ${Libxc_DIR}/lib/cmake/Libxc
    ${Libxc_DIR}/lib64/cmake/Libxc
  )
endif()

# Copy the results to the output variables and target.
# if find_package() above works, Libxc::xc would be present and
# below would be skipped.
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
