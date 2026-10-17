# CUDA adapter for the DDLA interface currently consumed by LibRI.
# DDLA_ROOT or CMAKE_PREFIX_PATH may name an installation prefix; the explicit
# DDLA_INCLUDE_DIR and DDLA_LIBRARY overrides remain supported.

find_path(DDLA_INCLUDE_DIR
  NAMES ddla/ddla_connector.h
  HINTS ${DDLA_ROOT}
  PATH_SUFFIXES include)
find_library(DDLA_LIBRARY
  NAMES ddla
  HINTS ${DDLA_ROOT}
  PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DDLA
  REQUIRED_VARS DDLA_INCLUDE_DIR DDLA_LIBRARY)

if(DDLA_FOUND AND NOT TARGET DDLA::DDLA)
  add_library(DDLA::DDLA UNKNOWN IMPORTED)
  set_target_properties(DDLA::DDLA PROPERTIES
    IMPORTED_LOCATION "${DDLA_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${DDLA_INCLUDE_DIR}"
    INTERFACE_COMPILE_DEFINITIONS DDLA_USE_CUDA
    INTERFACE_COMPILE_FEATURES cxx_std_17
    INTERFACE_LINK_LIBRARIES
      "MPI::MPI_CXX;CUDA::cudart;CUDA::cublas;CUDA::curand;CUDA::cusolver")
endif()

mark_as_advanced(DDLA_INCLUDE_DIR DDLA_LIBRARY)
