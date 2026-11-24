find_path(libnpy_SOURCE_DIR npy.hpp HINTS ${libnpy_INCLUDE_DIR})

if(NOT libnpy_SOURCE_DIR)
  include(FetchContent)
  FetchContent_Declare(
    libnpy
    GIT_REPOSITORY https://github.com/llohse/libnpy.git
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE)
  FetchContent_MakeAvailable(libnpy)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libnpy DEFAULT_MSG libnpy_SOURCE_DIR)

if(Libnpy_FOUND)
  if(NOT TARGET libnpy::npy)
    add_library(libnpy::npy UNKNOWN IMPORTED)
    set_target_properties(
      libnpy::npy
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${libnpy_SOURCE_DIR}/include")
  endif()
endif()
