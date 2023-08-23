add_library(
    operator_ks_lcao
    OBJECT
    ekinetic_lcao.cpp
    nonlocal_lcao.cpp
    op_exx_lcao.cpp
    op_dftu_lcao.cpp
    meta_lcao.cpp
    veff_lcao.cpp
    overlap_lcao.cpp
    deepks_lcao.cpp
    lambda_lcao.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(operator_ks_lcao)
endif()

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
