AddTest(
  TARGET MODULE_LCAO_extrap_test
  LIBS parameter ${math_libs} psi base device
  SOURCES
    wf_extrap_test.cpp
    ../wf_history_lcao.cpp
    ../wf_orthonormalize_lcao.cpp
    ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
)
