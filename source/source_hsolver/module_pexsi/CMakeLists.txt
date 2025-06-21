add_library(pexsi OBJECT dist_bcd_matrix.cpp dist_ccs_matrix.cpp dist_matrix_transformer.cpp pexsi_solver.cpp simple_pexsi.cpp)

if(ENABLE_COVERAGE)
  add_coverage(pexsi)
endif()
