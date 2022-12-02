list(APPEND objects
    fft.cpp
    pw_basis.cpp
    pw_basis_k.cpp
    pw_distributeg.cpp
    pw_distributeg_method1.cpp
    pw_distributeg_method2.cpp
    pw_distributer.cpp
    pw_init.cpp
    pw_transform.cpp
    pw_transform_k.cpp
    src/pw_multi_device.cpp
)

if (USE_CUDA)
  list(APPEND objects src/cuda/pw_multi_device.cu)
elseif(USE_ROCM)
  list(APPEND objects src/rocm/pw_multi_device.cu)
endif()

add_library(
    planewave
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(planewave)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
