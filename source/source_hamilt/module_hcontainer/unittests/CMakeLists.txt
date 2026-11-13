AddTest(
  TARGET MODULE_LCAO_hcontainer_transfer_test
  LIBS parameter psi base device
  SOURCES test_transfer.cpp
          ../func_transfer.cpp
          ../base_matrix.cpp
          ../hcontainer.cpp
          ../atom_pair.cpp
          ../transfer.cpp
          ../../../source_basis/module_ao/parallel_orbitals.cpp
          ../test/tmp_mocks.cpp
)
