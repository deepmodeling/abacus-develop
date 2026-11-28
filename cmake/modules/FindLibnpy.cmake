find_path(
  libnpy_INCLUDE_DIR npy.hpp
  HINTS ${libnpy_SOURCE_DIR}
  PATH_SUFFIXES include libnpy include/libnpy)
if(NOT libnpy_INCLUDE_DIR)
  include(FetchContent)
  FetchContent_Declare(
    libnpy
    GIT_REPOSITORY https://github.com/llohse/libnpy.git
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE)
  FetchContent_MakeAvailable(libnpy)
  find_path(
    libnpy_INCLUDE_DIR npy.hpp
    PATHS "${libnpy_SOURCE_DIR}/include"
    NO_DEFAULT_PATH)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libnpy DEFAULT_MSG libnpy_INCLUDE_DIR)

if(Libnpy_FOUND)
  if(NOT TARGET libnpy::npy)
    add_library(libnpy::npy INTERFACE IMPORTED)
    set_target_properties(libnpy::npy PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                 "${libnpy_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(libnpy_INCLUDE_DIR)
