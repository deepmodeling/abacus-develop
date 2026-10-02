abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(_OPENMP)

if (ENABLE_MPI)

AddTest(
  TARGET MODULE_DOS_cal_dos
  LIBS parameter base device
  SOURCES test_cal_dos.cpp ../cal_dos.cpp
)

AddTest(
  TARGET MODULE_DOS_write_dos_pw
  LIBS parameter base device symmetry
  SOURCES test_write_dos_pw.cpp ../write_dos_pw.cpp ../cal_dos.cpp
    ../../../source_cell/parallel_kpoints.cpp ../../../source_cell/klist.cpp
    ../../../source_cell/klist_io.cpp ../../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_DOS_write_pdos_text
  LIBS parameter base device
  SOURCES test_write_pdos_text.cpp ../write_pdos_text.cpp
    ../../../source_cell/atom_spec.cpp
    ../../../source_cell/atom_pseudo.cpp
    ../../../source_cell/pseudo.cpp
    ../../../source_cell/magnetism.cpp
    ../../../source_cell/test/support/mock_unitcell.cpp
)

endif()
