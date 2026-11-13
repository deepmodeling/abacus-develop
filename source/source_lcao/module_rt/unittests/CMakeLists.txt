AddTest(
  TARGET MODULE_LCAO_tddft_propagator_test
  LIBS parameter base device tddft_test_lib
  SOURCES test_propagator.cpp
          ../propagator.cpp
          ../propagator_cn2.cpp
          ../propagator_taylor.cpp
          ../propagator_etrs.cpp
          ../../../source_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET MODULE_LCAO_tddft_middle_hamilt_test
  LIBS parameter base device tddft_test_lib
  SOURCES test_middle_hamilt.cpp
          ../middle_hamilt.cpp
          ../../../source_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET MODULE_LCAO_tddft_upsi_test
  LIBS parameter base device tddft_test_lib
  SOURCES test_upsi.cpp
          ../upsi.cpp
          ../../../source_basis/module_ao/parallel_orbitals.cpp
)
