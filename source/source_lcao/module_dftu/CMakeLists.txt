list(APPEND objects
      dftu_nao.cpp
      dftu_nao_fs_k.cpp
      dftu_nao_folding.cpp
      dftu_nao_pots.cpp
      dftu_nao_occ.cpp
      dftu_nao_energy.cpp
      dftu_hamilt.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
        dftu_nao_op.cpp
        dftu_nao_fs_r.cpp
        dftu_nao_for_r.cpp
        dftu_nao_str_r.cpp
        dftu_nao_op_legacy.cpp
  )
endif()

add_library(
    dftu
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(dftu)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
