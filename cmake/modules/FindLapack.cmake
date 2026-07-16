# In compatibility to builtin FindLAPACK.cmake before v3.5.4
if(DEFINED LAPACK_DIR)
  string(APPEND CMAKE_PREFIX_PATH ";${LAPACK_DIR}")
endif()
if(DEFINED LAPACK_LIBRARY)
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
endif()

# find_package(Blas) must resolve to our cmake/FindBlas.cmake wrapper, so
# leave CMAKE_MODULE_PATH intact for it.
find_package(Blas REQUIRED)

# Delegate to CMake's builtin FindLAPACK module. As with FindBlas, the names
# "FindLapack.cmake" and builtin "FindLAPACK.cmake" collide on
# case-insensitive filesystems, so drop our module directory from
# CMAKE_MODULE_PATH around the call to avoid infinite recursion.
set(_abacus_lapack_saved_module_path "${CMAKE_MODULE_PATH}")
list(REMOVE_ITEM CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_package(LAPACK REQUIRED)
set(CMAKE_MODULE_PATH "${_abacus_lapack_saved_module_path}")

if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK UNKNOWN IMPORTED)
    set_target_properties(LAPACK::LAPACK PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${LAPACK_LIBRARIES}")
endif()
