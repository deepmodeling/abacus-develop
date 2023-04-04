add_library(
    psi
    OBJECT
    psi.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(psi)
endif()

if (BUILD_TESTING)
  add_subdirectory(kernels/test)
  add_subdirectory(test)
endif()