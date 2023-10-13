list(APPEND objects
    esolver.cpp
    esolver_ks.cpp
    esolver_fp.cpp
    esolver_ks_pw.cpp
    esolver_ks_pw_tool.cpp
    esolver_sdft_pw.cpp
    esolver_sdft_pw_tool.cpp
    esolver_lj.cpp
    esolver_dp.cpp
    esolver_of.cpp
)
if(ENABLE_LCAO)
  list(APPEND objects
      esolver_ks_lcao.cpp
      esolver_ks_lcao_elec.cpp
      esolver_ks_lcao_tddft.cpp
      esolver_ks_lcao_tmpfunc.cpp
  )
endif()

add_library(
    esolver
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(esolver)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

