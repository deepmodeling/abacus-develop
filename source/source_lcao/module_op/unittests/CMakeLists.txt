addtest(
  TARGET MODULE_LCAO_angmom_op_mat_test
  LIBS parameter base device neighbor hcontainer
  SOURCES test_angmom_op_mat.cpp
          ../angmom_op_mat.cpp
          ../../../source_io/module_hs/hs_sparse_io.cpp
          ../../../source_io/module_hs/lat_r_csr.cpp
          ../../../source_cell/ucell_io.cpp
          ../../../source_basis/module_ao/parallel_orbitals.cpp
          ../../../source_basis/module_ao/orb_atomic_lm.cpp
          ../../../source_basis/module_ao/orb_atomic.cpp
          ../../../source_basis/module_nao/radial_set.cpp
          ../../../source_basis/module_nao/numerical_radial.cpp
          ../../../source_basis/module_nao/beta_radials.cpp
          ../../../source_basis/module_nao/hydrogen_radials.cpp
          ../../../source_basis/module_nao/pswfc_radials.cpp
          ../../../source_basis/module_nao/atomic_radials.cpp
          ../../../source_basis/module_nao/smoothgen.cpp
          ../../../source_basis/module_nao/sphbes_radials.cpp
          ../../../source_basis/module_nao/two_center_integrator.cpp
          ../../../source_basis/module_nao/two_center_table.cpp
          ../../../source_basis/module_nao/real_gaunt_table.cpp
          ../../../source_basis/module_nao/radial_collection.cpp)

addtest(
  TARGET MODULE_LCAO_pos_op_csr_test
  LIBS parameter base device hcontainer
  SOURCES test_pos_op_csr.cpp
          ../pos_op_csr.cpp
          ../../../source_io/module_hs/lat_r_csr.cpp
          ../../../source_cell/ucell_io.cpp
          ../../../source_basis/module_ao/parallel_orbitals.cpp)
