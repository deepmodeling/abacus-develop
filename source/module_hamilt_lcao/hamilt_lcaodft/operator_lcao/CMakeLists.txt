add_library(
    operator_ks_lcao
    OBJECT
    op_exx_lcao.cpp
    op_dftu_lcao.cpp
    meta_lcao.cpp
    veff_lcao.cpp
    deepks_lcao.cpp
    overlap_new.cpp
    ekinetic_new.cpp
    nonlocal_new.cpp
    td_ekinetic_lcao.cpp
    td_nonlocal_lcao.cpp
    sc_lambda_lcao.cpp
    dftu_new.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(operator_ks_lcao)
endif()

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
