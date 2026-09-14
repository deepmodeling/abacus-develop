add_library(
    operator_ks_lcao
    OBJECT
    op_exx_lcao.cpp
    meta_lcao.cpp
    veff_lcao.cpp
    veff_dh.cpp
    deepks_lcao.cpp
    overlap.cpp
    overlap_fs.cpp
    ekinetic.cpp
    ekinetic_fs.cpp
    ekinetic_dh.cpp
    nonlocal.cpp
    nonlocal_dh.cpp
    nonlocal_fs.cpp
    td_ekinetic_lcao.cpp
    td_nonlocal_lcao.cpp
    td_pot_hybrid.cpp
    td_pot_hybrid_fs.cpp
    dspin_lcao.cpp
    dspin_fs.cpp
    operator_fs_utils.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(operator_ks_lcao)
endif()

IF (BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
