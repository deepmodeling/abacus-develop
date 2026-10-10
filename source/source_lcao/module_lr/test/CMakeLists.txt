if(ENABLE_LIBXC)
  AddTest(
    TARGET MODULE_LR_xc_kxc
    LIBS Libxc::xc
    SOURCES test_xc_kxc.cpp ../potentials/xc_kxc.cpp
  )
endif()

AddTest(
  TARGET MODULE_LR_exx_projection
  LIBS base parameter ${math_libs} container device
  SOURCES test_exx_proj.cpp ../exx_proj.cpp
)

AddTest(
  TARGET MODULE_LR_grad_degen
  LIBS base parameter ${math_libs} container device
  SOURCES test_grad_degen.cpp test_grad_jt.cpp ../grad_degen.cpp ../grad_jt.cpp
)

if(ENABLE_MPI)
  if(TARGET ELPA::ELPA)
    AddTest(
      TARGET MODULE_LR_zeqlin_solv
      LIBS base parameter ${math_libs} container device psi ELPA::ELPA MPI::MPI_CXX
      SOURCES test_zeqlin.cpp ../zeqlin_solv.cpp ../utils/lr_util.cpp
    )
  else()
    AddTest(
      TARGET MODULE_LR_zeqlin_solv
      LIBS base parameter ${math_libs} container device psi MPI::MPI_CXX
      SOURCES test_zeqlin.cpp ../zeqlin_solv.cpp ../utils/lr_util.cpp
    )
  endif()
endif()

AddTest(
  TARGET MODULE_LR_gradient_amplitudes
  LIBS base parameter ${math_libs} container device
  SOURCES test_lr_amp.cpp
)

# AO construction shares the existing two-center mock fixture with operator tests.
set(root_overlap_test_sources
  ../state_ovlp.cpp
  ../../module_operator_lcao/ovlp_block.cpp
  ../../module_operator_lcao/test/tmp_mocks.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/orb_atomic_lm.cpp
  ../../../source_hamilt/operator.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp
  ../../../source_hamilt/module_hcontainer/hcontainer.cpp
  ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_hamilt/module_hcontainer/func_folding.cpp
)

AddTest(
  TARGET MODULE_LR_state_ovlp
  LIBS base parameter ${math_libs} container device psi
  SOURCES test_state_ovlp.cpp ${root_overlap_test_sources}
)

AddTest(
  TARGET MODULE_LR_state_group
  LIBS base parameter ${math_libs} container device
  SOURCES test_state_grp.cpp ../state_group.cpp
)
