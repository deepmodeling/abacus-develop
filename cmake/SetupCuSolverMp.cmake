# =============================================================================
# Configure cuSOLVERMp dependencies and linking for ABACUS
# =============================================================================

include_guard(GLOBAL)

function(abacus_setup_cusolvermp target_name)
  add_compile_definitions(__CUSOLVERMP)

  # Find cuSOLVERMp first, then decide communicator backend.
  find_library(CUSOLVERMP_LIBRARY NAMES cusolverMp
      HINTS ${CAL_CUSOLVERMP_PATH} ${NVHPC_ROOT_DIR}
      PATH_SUFFIXES lib lib64 math_libs/lib math_libs/lib64)

  find_path(CUSOLVERMP_INCLUDE_DIR NAMES cusolverMp.h
      HINTS ${CAL_CUSOLVERMP_PATH} ${NVHPC_ROOT_DIR}
      PATH_SUFFIXES include math_libs/include)

  if(NOT CUSOLVERMP_LIBRARY OR NOT CUSOLVERMP_INCLUDE_DIR)
    message(FATAL_ERROR
      "cuSOLVERMp not found. Set CUSOLVERMP_PATH or NVHPC_ROOT_DIR."
    )
  endif()

  message(STATUS "Found cuSOLVERMp: ${CUSOLVERMP_LIBRARY}")

  set(CUSOLVERMP_VERSION_STR "")
  set(CUSOLVERMP_VERSION_HEADER "${CUSOLVERMP_INCLUDE_DIR}/cusolverMp.h")
  if(EXISTS "${CUSOLVERMP_VERSION_HEADER}")
    file(STRINGS "${CUSOLVERMP_VERSION_HEADER}" CUSOLVERMP_MAJOR_LINE
         REGEX "^#define[ \t]+CUSOLVERMP_VER_MAJOR[ \t]+[0-9]+")
    file(STRINGS "${CUSOLVERMP_VERSION_HEADER}" CUSOLVERMP_MINOR_LINE
         REGEX "^#define[ \t]+CUSOLVERMP_VER_MINOR[ \t]+[0-9]+")
    file(STRINGS "${CUSOLVERMP_VERSION_HEADER}" CUSOLVERMP_PATCH_LINE
         REGEX "^#define[ \t]+CUSOLVERMP_VER_PATCH[ \t]+[0-9]+")
    string(REGEX MATCH "([0-9]+)" CUSOLVERMP_VER_MAJOR "${CUSOLVERMP_MAJOR_LINE}")
    string(REGEX MATCH "([0-9]+)" CUSOLVERMP_VER_MINOR "${CUSOLVERMP_MINOR_LINE}")
    string(REGEX MATCH "([0-9]+)" CUSOLVERMP_VER_PATCH "${CUSOLVERMP_PATCH_LINE}")
    if(NOT CUSOLVERMP_VER_MAJOR STREQUAL ""
       AND NOT CUSOLVERMP_VER_MINOR STREQUAL ""
       AND NOT CUSOLVERMP_VER_PATCH STREQUAL "")
      set(CUSOLVERMP_VERSION_STR
          "${CUSOLVERMP_VER_MAJOR}.${CUSOLVERMP_VER_MINOR}.${CUSOLVERMP_VER_PATCH}")
    endif()
  endif()

  # Check minimum version requirement (>= 0.4.0)
  if(CUSOLVERMP_VERSION_STR AND CUSOLVERMP_VERSION_STR VERSION_LESS "0.4.0")
    message(FATAL_ERROR
      "cuSOLVERMp version ${CUSOLVERMP_VERSION_STR} is too old. "
      "ABACUS requires cuSOLVERMp >= 0.4.0 (NVIDIA HPC SDK >= 23.5). "
      "Please upgrade your NVIDIA HPC SDK installation."
    )
  endif()

  # Auto-select communicator backend by cuSOLVERMp version.
  # cuSOLVERMp < 0.7.0 -> CAL, otherwise -> NCCL.
  set(_use_cal OFF)
  if(CUSOLVERMP_VERSION_STR AND CUSOLVERMP_VERSION_STR VERSION_LESS "0.7.0")
    set(_use_cal ON)
    message(STATUS
      "Detected cuSOLVERMp ${CUSOLVERMP_VERSION_STR} (< 0.7.0). Using CAL backend.")
  elseif(CUSOLVERMP_VERSION_STR)
    message(STATUS
      "Detected cuSOLVERMp ${CUSOLVERMP_VERSION_STR} (>= 0.7.0). Using NCCL backend.")
  elseif(NOT CUSOLVERMP_VERSION_STR)
    message(WARNING
      "Unable to detect cuSOLVERMp version from header. Using NCCL backend by default.")
  endif()

  # Raise the variable to the caller's scope
  set(_use_cal ${_use_cal} PARENT_SCOPE)

  # Backend selection:
  # - _use_cal=ON  -> cal communicator backend
  # - _use_cal=OFF -> NCCL communicator backend
  if(_use_cal)
    add_compile_definitions(__USE_CAL)

    find_library(CAL_LIBRARY NAMES cal
        HINTS ${CAL_CUSOLVERMP_PATH} ${NVHPC_ROOT_DIR}
        PATH_SUFFIXES lib lib64 math_libs/lib64)
    find_path(CAL_INCLUDE_DIR NAMES cal.h
        HINTS ${CAL_CUSOLVERMP_PATH} ${NVHPC_ROOT_DIR}
        PATH_SUFFIXES include math_libs/include)

    if(NOT CAL_LIBRARY OR NOT CAL_INCLUDE_DIR)
      message(FATAL_ERROR "CAL not found. Set CAL_PATH or NVHPC_ROOT_DIR.")
    endif()

    message(STATUS "Found CAL: ${CAL_LIBRARY}")
    if(NOT TARGET CAL::CAL)
      add_library(CAL::CAL IMPORTED INTERFACE)
      set_target_properties(CAL::CAL PROPERTIES
          INTERFACE_LINK_LIBRARIES "${CAL_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${CAL_INCLUDE_DIR}")
    endif()
  else()

    find_library(NCCL_LIBRARY NAMES nccl
        HINTS ${NCCL_PATH} ${NVHPC_ROOT_DIR}
        PATH_SUFFIXES lib lib64 comm_libs/nccl/lib)
    find_path(NCCL_INCLUDE_DIR NAMES nccl.h
        HINTS ${NCCL_PATH} ${NVHPC_ROOT_DIR}
        PATHS ${CUDA_TOOLKIT_ROOT_DIR}
        PATH_SUFFIXES include comm_libs/nccl/include)

    if(NOT NCCL_LIBRARY OR NOT NCCL_INCLUDE_DIR)
      message(FATAL_ERROR "NCCL not found. Set NCCL_PATH or NVHPC_ROOT_DIR.")
    endif()

    message(STATUS "Found NCCL: ${NCCL_LIBRARY}")
    if(NOT TARGET NCCL::NCCL)
      add_library(NCCL::NCCL IMPORTED INTERFACE)
      set_target_properties(NCCL::NCCL PROPERTIES
          INTERFACE_LINK_LIBRARIES "${NCCL_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${NCCL_INCLUDE_DIR}")
    endif()
  endif()

  # Create cusolverMp::cusolverMp imported target
  if(NOT TARGET cusolverMp::cusolverMp)
    add_library(cusolverMp::cusolverMp IMPORTED INTERFACE)
    set_target_properties(cusolverMp::cusolverMp PROPERTIES
        INTERFACE_LINK_LIBRARIES "${CUSOLVERMP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CUSOLVERMP_INCLUDE_DIR}")
  endif()

  # === Link libraries ===
  if(_use_cal)
    target_link_libraries(${target_name}
        CAL::CAL
        cusolverMp::cusolverMp)
  else()
    target_link_libraries(${target_name}
        NCCL::NCCL
        cusolverMp::cusolverMp)
  endif()
endfunction()
