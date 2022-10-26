list(APPEND psi_srcs
    psi.cpp
    src/device.cpp
    src/memory_psi.cpp
)

if (USE_CUDA) 
    list(APPEND psi_srcs src/cuda/memory.cu)
elseif(USE_ROCM)
    list(APPEND psi_srcs src/rocm/memory.hip)
endif()

add_library(
    psi
    OBJECT
    ${psi_srcs}
)

if(ENABLE_COVERAGE)
  add_coverage(psi)
endif()

if (BUILD_TESTING)
  add_subdirectory(test)
endif()