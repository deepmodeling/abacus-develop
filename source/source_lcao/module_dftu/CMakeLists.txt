# The operators below are LCAO-specific. dftu_dc_check.cpp is compiled by
# io_input so INPUT validation remains available in builds without LCAO.
# The Plus_U_Base class lives in source_pw/module_pwdft/dftu_base.* and is
# compiled unconditionally there.
if(ENABLE_LCAO)
  list(
    APPEND
    objects
    dftu_nao_fs_k.cpp
    dftu_nao_folding.cpp
    dftu_nao_pots.cpp
    dftu_full_dc.cpp
    dftu_nao_occ.cpp
    dftu_nao_energy.cpp
    dftu_nao_op.cpp
    dftu_nao_adj.cpp
    dftu_nao_fs_r.cpp
    dftu_nao_fs_reduce.cpp
    dftu_nao_for_r.cpp
    dftu_nao_str_r.cpp
    dftu_nao_op_legacy.cpp)
  add_library(dftu OBJECT ${objects})

  if(ENABLE_COVERAGE)
    add_coverage(dftu)
  endif()

  if(BUILD_TESTING)
    add_subdirectory(unittests)
  endif()
endif()
