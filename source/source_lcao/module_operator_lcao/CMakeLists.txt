add_library(
    operator_ks_lcao
    OBJECT
    op_exx_lcao.cpp
    op_dftu_lcao.cpp
    meta_lcao.cpp
    veff_lcao.cpp
    deepks_lcao.cpp
    overlap.cpp
    ekinetic.cpp
    nonlocal.cpp
    td_ekinetic_lcao.cpp
    td_nonlocal_lcao.cpp
    td_pot_hybrid.cpp
    dspin_lcao.cpp
    dftu_lcao.cpp
    operator_force_stress_utils.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(operator_ks_lcao)
endif()

IF (BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
