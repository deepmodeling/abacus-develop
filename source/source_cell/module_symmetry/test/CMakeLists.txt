abacus_disable_feature_definitions(__LCAO)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
AddTest(
  TARGET MODULE_CELL_SYMMETRY_analysis
  LIBS base device symmetry
  SOURCES symm_test.cpp symm_test_analysis.cpp
)
AddTest(
  TARGET MODULE_CELL_SYMMETRY_symtrz
  LIBS base device symmetry
  SOURCES symm_test.cpp symm_test_symtrz.cpp
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
    ${ABACUS_SOURCE_DIR}/source_estate/module_dm/density_matrix.cpp
    ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
    ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
    ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
    ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
)