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
  ../../../source_cell/klist.cpp ../../../source_cell/parallel_kpoints.cpp ../../../source_cell/reciprocal_grid.cpp
)