find_package(ZLIB REQUIRED)

find_path(
  CNPY_INCLUDE_DIR cnpy.h
  HINTS ${CNPY_DIR} ${cnpy_SOURCE_DIR}
  PATH_SUFFIXES include)
find_library(
  CNPY_LIBRARY
  NAMES cnpy
  HINTS ${CNPY_DIR} "${CNPY_INCLUDE_DIR}/.."
  PATH_SUFFIXES lib lib64)

if(NOT CNPY_INCLUDE_DIR OR NOT CNPY_LIBRARY)
  include(FetchContent)
  FetchContent_Declare(
    cnpy
    GIT_REPOSITORY https://github.com/rogersce/cnpy.git
    GIT_PROGRESS TRUE)
  FetchContent_MakeAvailable(cnpy)
  set(CNPY_INCLUDE_DIR "${cnpy_SOURCE_DIR}")
  set(CNPY_LIBRARY cnpy)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CNPY DEFAULT_MSG CNPY_INCLUDE_DIR
                                  CNPY_LIBRARY)

if(CNPY_FOUND AND NOT TARGET CNPY::cnpy)
  add_library(CNPY::cnpy INTERFACE IMPORTED)
  set_target_properties(
    CNPY::cnpy PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CNPY_INCLUDE_DIR}"
                          INTERFACE_LINK_LIBRARIES "${CNPY_LIBRARY};ZLIB::ZLIB")
endif()

mark_as_advanced(CNPY_INCLUDE_DIR CNPY_LIBRARY)
