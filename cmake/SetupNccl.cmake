include_guard(GLOBAL)

include(CheckIncludeFileCXX)

function(abacus_setup_nccl target_name)
  find_library(NCCL_LIBRARY NAMES nccl
      HINTS ${NCCL_PATH} ${NVHPC_ROOT_DIR}
      PATH_SUFFIXES lib lib64 comm_libs/nccl/lib)
  find_path(NCCL_INCLUDE_DIR NAMES nccl.h
      HINTS ${NCCL_PATH} ${NVHPC_ROOT_DIR}
      PATHS ${CUDAToolkit_ROOT}
      PATH_SUFFIXES include comm_libs/nccl/include)

  check_include_file_cxx("nccl.h" HAVE_NCCL_HEADER)

  if(NOT NCCL_LIBRARY)
    set(NCCL_LIBRARY nccl)
  endif()

  if(NOT NCCL_INCLUDE_DIR AND NOT HAVE_NCCL_HEADER)
    message(FATAL_ERROR
      "NCCL not found. Set NCCL_PATH or NVHPC_ROOT_DIR.")
  endif()

  if(NCCL_INCLUDE_DIR)
    message(STATUS "Found NCCL for parallel_device: ${NCCL_LIBRARY}")
  else()
    message(STATUS "Using default compiler/linker search paths for NCCL: ${NCCL_LIBRARY}")
  endif()
  if(NOT TARGET NCCL::NCCL)
    add_library(NCCL::NCCL IMPORTED INTERFACE)
    if(NCCL_INCLUDE_DIR)
      set_target_properties(NCCL::NCCL PROPERTIES
          INTERFACE_LINK_LIBRARIES "${NCCL_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${NCCL_INCLUDE_DIR}")
    else()
      set_target_properties(NCCL::NCCL PROPERTIES
          INTERFACE_LINK_LIBRARIES "${NCCL_LIBRARY}")
    endif()
  endif()

  if(NCCL_INCLUDE_DIR)
    # `parallel_device.cpp` is compiled inside the later `base` OBJECT library,
    # so the header path must also be visible to targets created in subdirs.
    include_directories(${NCCL_INCLUDE_DIR})
    target_include_directories(${target_name} PRIVATE ${NCCL_INCLUDE_DIR})
  endif()
  target_link_libraries(${target_name} NCCL::NCCL)
endfunction()
