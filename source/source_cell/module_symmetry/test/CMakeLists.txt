abacus_disable_feature_definitions(__LCAO)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
AddTest(
  TARGET MODULE_CELL_SYMMETRY_analysis
  LIBS base device symmetry
  SOURCES symmetry_test.cpp symmetry_test_analysis.cpp
)
AddTest(
  TARGET MODULE_CELL_SYMMETRY_symtrz
  LIBS base device symmetry
  SOURCES symmetry_test.cpp symmetry_test_symtrz.cpp
)
AddTest(
  TARGET MODULE_CELL_SYMMETRY_rotation_spin
  LIBS parameter base ${math_libs} device
  SOURCES symm_rot_spin_test.cpp ../symm_rot_spin.cpp
)
AddTest(
  TARGET MODULE_CELL_SYMMETRY_rho_soc
  LIBS parameter base ${math_libs} device symmetry
  SOURCES symm_rho_soc_test.cpp
)