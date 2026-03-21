# =============================================================================
# Configure cuBLASMp dependencies and linking for ABACUS
# =============================================================================

include_guard(GLOBAL)

function(abacus_setup_cublasmp target_name)
  add_compile_definitions(__CUBLASMP)

  # 1. Search for cuBLASMp library and header files
  # libcublasmp.so
  find_library(CUBLASMP_LIBRARY NAMES cublasmp
      HINTS ${CUBLASMP_PATH} ${NVHPC_ROOT_DIR}
      PATH_SUFFIXES lib lib64 math_libs/lib math_libs/lib64)

  # cublasmp.h
  find_path(CUBLASMP_INCLUDE_DIR NAMES cublasmp.h
      HINTS ${CUBLASMP_PATH} ${NVHPC_ROOT_DIR}
      PATH_SUFFIXES include math_libs/include)

  if(NOT CUBLASMP_LIBRARY OR NOT CUBLASMP_INCLUDE_DIR)
    message(FATAL_ERROR
      "cuBLASMp not found. Please ensure CUBLASMP_PATH is set correctly."
    )
  endif()

  message(STATUS "Found cuBLASMp: ${CUBLASMP_LIBRARY}")

  # 2. Version validation by parsing header macros
  set(CUBLASMP_VERSION_STR "")
  set(CUBLASMP_VERSION_HEADER "${CUBLASMP_INCLUDE_DIR}/cublasmp.h")
  
  if(EXISTS "${CUBLASMP_VERSION_HEADER}")
    # Extract version lines using regular expressions from cublasmp.h
    file(STRINGS "${CUBLASMP_VERSION_HEADER}" CUBLASMP_MAJOR_LINE
         REGEX "^#define[ \t]+CUBLASMP_VER_MAJOR[ \t]+[0-9]+")
    file(STRINGS "${CUBLASMP_VERSION_HEADER}" CUBLASMP_MINOR_LINE
         REGEX "^#define[ \t]+CUBLASMP_VER_MINOR[ \t]+[0-9]+")
    file(STRINGS "${CUBLASMP_VERSION_HEADER}" CUBLASMP_PATCH_LINE
         REGEX "^#define[ \t]+CUBLASMP_VER_PATCH[ \t]+[0-9]+")
    
    # Extract numeric values from the matched strings
    string(REGEX MATCH "([0-9]+)" CUBLASMP_VER_MAJOR "${CUBLASMP_MAJOR_LINE}")
    string(REGEX MATCH "([0-9]+)" CUBLASMP_VER_MINOR "${CUBLASMP_MINOR_LINE}")
    string(REGEX MATCH "([0-9]+)" CUBLASMP_VER_PATCH "${CUBLASMP_PATCH_LINE}")
    
    if(NOT CUBLASMP_VER_MAJOR STREQUAL ""
       AND NOT CUBLASMP_VER_MINOR STREQUAL ""
       AND NOT CUBLASMP_VER_PATCH STREQUAL "")
      set(CUBLASMP_VERSION_STR
          "${CUBLASMP_VER_MAJOR}.${CUBLASMP_VER_MINOR}.${CUBLASMP_VER_PATCH}")
    endif()
  endif()

  message(STATUS "Detected cuBLASMp version: ${CUBLASMP_VERSION_STR}")

  # 3. Version constraint: ABACUS requires cuBLASMp >= 0.8.0
  if(CUBLASMP_VERSION_STR AND CUBLASMP_VERSION_STR VERSION_LESS "0.8.0")
    message(FATAL_ERROR
      "cuBLASMp version ${CUBLASMP_VERSION_STR} is too old. "
      "ABACUS requires cuBLASMp >= 0.8.0 for NCCL Symmetric Memory support."
    )
  elseif(NOT CUBLASMP_VERSION_STR)
    message(WARNING "Could not detect cuBLASMp version. Proceeding cautiously.")
  endif()

  # 4. Create cublasMp::cublasMp imported target
  if(NOT TARGET cublasMp::cublasMp)
    add_library(cublasMp::cublasMp IMPORTED INTERFACE)
    set_target_properties(cublasMp::cublasMp PROPERTIES
        INTERFACE_LINK_LIBRARIES "${CUBLASMP_LIBRARY};NCCL::NCCL"
        INTERFACE_INCLUDE_DIRECTORIES "${CUBLASMP_INCLUDE_DIR}")
  endif()

  # 5. Link the library to the target
  target_link_libraries(${target_name} cublasMp::cublasMp)

endfunction()
