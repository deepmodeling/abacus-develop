add_library(
    psi
    OBJECT
    psi.cpp
)

add_library(
    psi_overall_init
    OBJECT
    setup_psi_pw.cpp
    setup_psi.cpp
    psi_prepare.cpp
)

add_library(
  psi_initializer
  OBJECT
  psi_initializer.cpp
  psi_init_random.cpp
  psi_init_file.cpp
  psi_init_atomic.cpp
  psi_init_atomic_random.cpp
  psi_init_nao.cpp
  psi_init_nao_random.cpp
)


if(ENABLE_COVERAGE)
  add_coverage(psi)
  add_coverage(psi_initializer)
endif()

if (BUILD_TESTING)
  # add_subdirectory(kernels/test)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
