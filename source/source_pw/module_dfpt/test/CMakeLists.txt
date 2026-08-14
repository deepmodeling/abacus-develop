abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

AddTest(
  TARGET MODULE_DFPT_irrep_data_test
  LIBS parameter base device symmetry
  SOURCES dfpt_irrep_data_test.cpp
  	../dfpt_irrep_data.cpp
  	../dfpt_pw_data.cpp
  	../../../source_cell/qlist.cpp
  	../../../source_cell/reciprocal_grid.cpp
  	../../../source_psi/psi.cpp
)

AddTest(
  TARGET MODULE_DFPT_pw_run_test
  LIBS parameter base device symmetry
  SOURCES dfpt_pw_run_test.cpp
  	../dfpt_pw.cpp
  	../dfpt_pw_data.cpp
  	../dfpt_irrep_data.cpp
  	../dfpt_pert.cpp
  	../dfpt_stern.cpp
  	../dfpt_rho.cpp
  	../dfpt_phon.cpp
  	../dfpt_q0.cpp
  	../dfpt_metal.cpp
  	../../../source_cell/qlist.cpp
  	../../../source_cell/reciprocal_grid.cpp
  	../../../source_psi/psi.cpp
)
