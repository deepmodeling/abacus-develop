list(APPEND objects
    fft.cpp
    pw_basis.cpp
    pw_basis_k.cpp
    pw_basis_sup.cpp
    pw_distributeg.cpp
    pw_distributeg_method1.cpp
    pw_distributeg_method2.cpp
    pw_distributer.cpp
    pw_init.cpp
    pw_transform.cpp
    pw_transform_k.cpp
)

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
  add_subdirectory(test_serial)
  add_subdirectory(kernels/test)
endif()
