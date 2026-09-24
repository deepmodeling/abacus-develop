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
  TARGET MODULE_DFTU_pots
  LIBS base device
  SOURCES test_dftu_nao_pots.cpp
)

AddTest(
  TARGET MODULE_DFTU_energy
  LIBS base device
  SOURCES test_dftu_nao_energy.cpp
)

AddTest(
  TARGET MODULE_DFTU_op_legacy
  LIBS base device
  SOURCES test_dftu_nao_op_legacy.cpp
)

AddTest(
  TARGET MODULE_DFTU_for_r
  LIBS base device
  SOURCES test_dftu_nao_for_r.cpp
)

AddTest(
  TARGET MODULE_DFTU_str_r
  LIBS base device
  SOURCES test_dftu_nao_str_r.cpp
)

AddTest(
  TARGET MODULE_DFTU_folding
  LIBS base device
  SOURCES test_dftu_nao_folding.cpp
)

AddTest(
  TARGET MODULE_DFTU_fs_r
  LIBS base device
  SOURCES test_dftu_nao_fs_r.cpp ../dftu_nao_fs_reduce.cpp
)

AddTest(
  TARGET MODULE_DFTU_fs_reduce
  LIBS base device
  SOURCES test_dftu_nao_fs_reduce.cpp ../dftu_nao_fs_reduce.cpp
)

AddTest(
  TARGET MODULE_DFTU_fs_k
  LIBS base device
  SOURCES test_dftu_nao_fs_k.cpp
)

AddTest(
  TARGET MODULE_DFTU_fs_accum
  LIBS base device
  SOURCES test_dftu_nao_fs_accum.cpp ../../../source_basis/module_ao/parallel_orbitals.cpp
  # This test uses Parallel_2D::set_serial which is a real implementation
  # in parallel_2d.cpp; keep __MPI so the serial path links correctly.
  KEEP_FEATURE_DEFINITIONS __MPI
)

if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_DFTU_op
  LIBS parameter psi base device container symmetry
  SOURCES test_dftu_nao_op.cpp ../dftu_nao_op.cpp ../dftu_nao_adj.cpp ../dftu_nao_pots.cpp ../dftu_nao_fs_r.cpp ../dftu_nao_fs_reduce.cpp ../dftu_nao_for_r.cpp ../dftu_nao_str_r.cpp  ../../../source_estate/module_dm/density_matrix.cpp
  ../../../source_estate/module_dm/dmr_gamma.cpp
  ../../../source_estate/module_dm/dmr_init.cpp
  ../../../source_estate/module_dm/dm_setter.cpp
  ../../../source_estate/module_dm/dm_getter.cpp
  ../../../source_estate/module_dm/dm_tools.cpp
  ../../../source_estate/module_dm/dmr_k.cpp
  ../../../source_estate/module_dm/dmr_td.cpp
  ../../../source_estate/module_dm/dmr_full.cpp
  ../../../source_pw/module_pwdft/dftu_base.cpp
  ../../../source_pw/module_pwdft/dftu_base_io.cpp
  ../../../source_pw/module_pwdft/yukawa_screening.cpp
  ../../../source_estate/occ_matrix.cpp
  ../../../source_estate/occ_mixer.cpp
  ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp
  ../../../source_hamilt/module_hcontainer/hcontainer.cpp
  ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/orb_atomic_lm.cpp
  ../../module_operator_lcao/test/tmp_mocks.cpp
  ../../../source_hamilt/operator.cpp
  # This test calls MPI_Init in main() and drives ScaLAPACK distribution,
  # so it must keep __MPI even though this directory disables it.
  KEEP_FEATURE_DEFINITIONS __MPI
)

AddTest(
  TARGET MODULE_DFTU_ijr
  LIBS parameter psi base device container
  SOURCES test_dftu_nao_ijr.cpp
  ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp
  ../../../source_hamilt/module_hcontainer/hcontainer.cpp
  ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/orb_atomic_lm.cpp
  ../../module_operator_lcao/test/tmp_mocks.cpp
  ../../../source_hamilt/operator.cpp
  KEEP_FEATURE_DEFINITIONS __MPI
)
endif()

endif()
