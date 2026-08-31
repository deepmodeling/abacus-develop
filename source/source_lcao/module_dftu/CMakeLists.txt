list(APPEND objects
      dftu_lcao.cpp
      dftu_force.cpp
      dftu_yukawa.cpp
      dftu_folding.cpp
      dftu_lcao_pots.cpp
      dftu_lcao_occ.cpp
      dftu_lcao_energy.cpp
      dftu_hamilt.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
        dftu_lcao_op.cpp
        dftu_fs.cpp
        dftu_lcao_op_legacy.cpp
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
