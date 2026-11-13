abacus_disable_feature_definitions(__MPI)
# These tests exercise CPU code only and link no GPU kernels, so build them
# without the device instantiations (as source_hsolver/unittests does).
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

addtest(
  TARGET MODULE_PW_Sto_Tool_UTs
  LIBS parameter psi base device
  SOURCES ../sto_tool.cpp test_sto_tool.cpp)

addtest(
  TARGET MODULE_PW_Sto_Hamilt_UTs
  LIBS parameter psi base device planewave_serial symmetry
  SOURCES ../sto_hamilt_pw.cpp
          test_sto_hamilt_pw.cpp
          ../../../source_hamilt/operator.cpp
          ../../../source_cell/klist.cpp
          ../../../source_cell/klist_io.cpp
          ../../../source_cell/parallel_kpoints.cpp
          ../../../source_cell/reciprocal_grid.cpp)
