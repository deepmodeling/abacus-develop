add_library(
    psi
    OBJECT
    psi.cpp
)

add_library(
  psi_initializer
  OBJECT
  psi_initializer.cpp
  psi_initializer_random.cpp
  psi_initializer_atomic.cpp
  psi_initializer_atomic_random.cpp
  psi_initializer_nao.cpp
  psi_initializer_nao_random.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(psi)
  add_coverage(psi_initializer)
endif()

if (BUILD_TESTING)
  add_subdirectory(kernels/test)
  add_subdirectory(test)

  add_compile_definitions(PSI_INITIALIZER_TEST)
endif()