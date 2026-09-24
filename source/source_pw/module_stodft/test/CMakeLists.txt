abacus_disable_feature_definitions(__MPI)

AddTest(
  TARGET MODULE_PW_Sto_Tool_UTs
  LIBS parameter psi base device
  SOURCES ../sto_tool.cpp test_sto_tool.cpp
)

AddTest(
  TARGET MODULE_PW_Sto_Hamilt_UTs
  LIBS parameter psi base device planewave_serial symmetry
  SOURCES ../hamilt_sdft_pw.cpp test_hamilt_sto.cpp ../../../source_hamilt/operator.cpp
  ../../../source_cell/klist.cpp ../../../source_cell/klist_io.cpp ../../../source_cell/parallel_kpoints.cpp ../../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_PW_Sto_HSolver_UTs
  LIBS parameter psi device base container MPI::MPI_CXX
  SOURCES test_hsolver_sdft_pw.cpp ../hsolver_sdft_pw.cpp
  ../../../source_hsolver/hsolver_pw.cpp ../../../source_hsolver/diago_bpcg.cpp
  ../../../source_hsolver/diago_dav_subspace.cpp ../../../source_hsolver/diag_const_nums.cpp
  ../../../source_hsolver/diago_iter_assist.cpp ../../../source_hsolver/para_lin_tf.cpp
  ../../../source_estate/elecstate_tools.cpp ../../../source_estate/occupy.cpp
  ../../../source_base/module_fft/fft_bundle.cpp ../../../source_base/module_fft/fft_cpu.cpp
  # This test calls MPI_Init in main() and its mocks take MPI_Comm arguments
  # unconditionally, so it must keep __MPI even though this directory disables it.
  KEEP_FEATURE_DEFINITIONS __MPI
)