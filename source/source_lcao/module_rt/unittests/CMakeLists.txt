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
