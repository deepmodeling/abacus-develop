if(DEFINED BLAS_DIR)
    string(APPEND CMAKE_PREFIX_PATH ";${BLAS_DIR}")
endif()
if(DEFINED BLAS_LIBRARY)
    set(BLAS_LIBRARIES ${BLAS_LIBRARY})
endif()

# Delegate to CMake's builtin FindBLAS module. On case-insensitive
# filesystems (Windows, macOS) this file "FindBlas.cmake" and the builtin
# "FindBLAS.cmake" resolve to the same name, so a plain find_package(BLAS)
# recurses into this very file. Temporarily remove our module directory
# from CMAKE_MODULE_PATH so the builtin module is used instead. Harmless
# no-op on case-sensitive filesystems.
set(_abacus_blas_saved_module_path "${CMAKE_MODULE_PATH}")
list(REMOVE_ITEM CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_package(BLAS REQUIRED)
set(CMAKE_MODULE_PATH "${_abacus_blas_saved_module_path}")

if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS UNKNOWN IMPORTED)
    set_target_properties(BLAS::BLAS PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
	IMPORTED_LOCATION "${BLAS_LIBRARIES}")
endif()
