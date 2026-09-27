# Find an external GreenX minimax installation or source tree.
#
# Supported inputs:
#   GREENX_DIR=<prefix>  install prefix containing GreenXConfig.cmake or
#                        Fortran modules and libgx_minimax
#   GREENX_DIR=<package> GreenX install prefix or package directory
#
# The module exposes GreenX::MiniMax when the minimax API is available.

set(GreenX_FOUND FALSE)

if(TARGET GreenX::MiniMax)
  set(GreenX_FOUND TRUE)
elseif(TARGET greenX::LibGXMiniMax)
  add_library(GreenX::MiniMax ALIAS greenX::LibGXMiniMax)
  set(GreenX_FOUND TRUE)
elseif(TARGET LibGXMiniMax)
  add_library(GreenX::MiniMax ALIAS LibGXMiniMax)
  set(GreenX_FOUND TRUE)
endif()

# GreenX is an external dependency.  Prefer its installed CMake package so
# that ABACUS does not compile or vendor GreenX sources.  GREENX_DIR may point
# either to the install prefix or directly to the package directory.
if(NOT GreenX_FOUND)
  set(_greenx_config_paths)
  if(GREENX_DIR)
    list(APPEND _greenx_config_paths
      "${GREENX_DIR}"
      "${GREENX_DIR}/lib/cmake/greenX"
      "${GREENX_DIR}/lib64/cmake/greenX"
      "${GREENX_DIR}/share/greenX/cmake")
  endif()
  if(GREENX_DIR)
    find_package(greenX CONFIG QUIET
      PATHS ${_greenx_config_paths}
      NO_DEFAULT_PATH)
  else()
    find_package(greenX CONFIG QUIET)
  endif()
  if(TARGET greenX::LibGXMiniMax)
    add_library(GreenX::MiniMax ALIAS greenX::LibGXMiniMax)
    set(GreenX_FOUND TRUE)
  endif()
endif()

if(NOT GreenX_FOUND AND GREENX_DIR)
  find_path(GreenX_Fortran_MODULE_DIR
    NAMES gx_minimax.mod
    PATHS
      "${GREENX_DIR}/include"
      "${GREENX_DIR}/modules"
      "${GREENX_DIR}/lib"
    NO_DEFAULT_PATH)
  find_library(GreenX_MINIMAX_LIBRARY
    NAMES gx_minimax LibGXMiniMax
    PATHS
      "${GREENX_DIR}/lib"
      "${GREENX_DIR}/lib64"
    NO_DEFAULT_PATH)
  find_library(GreenX_COMMON_LIBRARY
    NAMES GXCommon gx_common
    PATHS
      "${GREENX_DIR}/lib"
      "${GREENX_DIR}/lib64"
    NO_DEFAULT_PATH)

  if(GreenX_Fortran_MODULE_DIR AND GreenX_MINIMAX_LIBRARY)
    add_library(GreenX::MiniMax UNKNOWN IMPORTED)
    set_target_properties(GreenX::MiniMax PROPERTIES
      IMPORTED_LOCATION "${GreenX_MINIMAX_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GreenX_Fortran_MODULE_DIR}")
    if(math_libs)
      set_property(TARGET GreenX::MiniMax APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES "${math_libs}")
    endif()
    if(GreenX_COMMON_LIBRARY)
      set_property(TARGET GreenX::MiniMax APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES "${GreenX_COMMON_LIBRARY}")
    endif()
    set(GreenX_FOUND TRUE)
  endif()
endif()

if(NOT GreenX_FOUND AND GreenX_FIND_REQUIRED)
  message(FATAL_ERROR
    "GreenX minimax was requested but no usable installation was found. "
    "Set GREENX_DIR to a GreenX installation prefix or CMake package directory "
    "containing gx_minimax.mod and libgx_minimax.")
endif()

mark_as_advanced(GreenX_Fortran_MODULE_DIR GreenX_MINIMAX_LIBRARY GreenX_COMMON_LIBRARY)
