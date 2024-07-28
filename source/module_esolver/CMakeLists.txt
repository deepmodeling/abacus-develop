list(APPEND objects
    esolver.cpp
    esolver_ks.cpp
    esolver_fp.cpp
    esolver_ks_pw.cpp
    esolver_ks_lcaopw.cpp
    esolver_sdft_pw.cpp
    esolver_lj.cpp
    esolver_dp.cpp
    esolver_of.cpp
    esolver_of_interface.cpp
    esolver_of_tool.cpp
    print_funcs.cpp
)
if(ENABLE_LCAO)
  list(APPEND objects
      esolver_ks_lcao.cpp
      esolver_ks_lcao_elec.cpp
      esolver_ks_lcao_tddft.cpp
      dpks_cal_e_delta_band.cpp
      io_npz.cpp
      set_matrix_grid.cpp
      dftu_cal_occup_m.cpp
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
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()

