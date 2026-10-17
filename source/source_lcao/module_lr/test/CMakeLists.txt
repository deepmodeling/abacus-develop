if(ENABLE_MPI)
  set(test_zeqlin_solv_main ../../../source_base/test/mpi_test_main.cpp)
endif()
if(ENABLE_LIBXC)
  AddTest(
    TARGET MODULE_LR_pot_hxc_lrtd
    LIBS base parameter device container psi Libxc::xc
    SOURCES test_pot_hxc_lrtd.cpp ../potentials/pot_hxc_lrtd.cpp ../kernel_pots.cpp ../potentials/pot_grad_xc.cpp
      ../potentials/xc_kernel.cpp ../potentials/xc_kxc.cpp ../utils/lr_util.cpp
      ../../../source_hamilt/module_xc/xc_functional.cpp
      ../../../source_hamilt/module_xc/libxc_setup.cpp
  )

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
  AddTest(
    TARGET MODULE_LR_zeqlin_solv
    LIBS base parameter ${math_libs} container device psi MPI::MPI_CXX
    SOURCES test_zeqlin_solv.cpp ${test_zeqlin_solv_main} ../zeq_cg.cpp ../zeqlin_solv.cpp ../utils/lr_util.cpp
  )
  add_test(
    NAME MODULE_LR_zeqlin_failure_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:MODULE_LR_zeqlin_solv> ${MPIEXEC_POSTFLAGS}
      --gtest_filter=ZeqLinearSolverTest.*FailsCollectively*
  )
  set_tests_properties(MODULE_LR_zeqlin_failure_mpi PROPERTIES
    PROCESSORS 4 TIMEOUT 60
    ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
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
  SOURCES test_state_ovlp.cpp state_ovlp_test.h ${root_overlap_test_sources}
)

AddTest(
  TARGET MODULE_LR_state_group
  LIBS base parameter ${math_libs} container device
  SOURCES test_state_group.cpp ../state_group.cpp
)

set(state_track_test_sources test_state_track.cpp ../state_track.cpp ../state_group.cpp ../grad_degen.cpp
                             ${root_overlap_test_sources})
if(ENABLE_MPI)
  list(APPEND state_track_test_sources ../../../source_base/test/mpi_test_main.cpp)
endif()
AddTest(
  TARGET MODULE_LR_state_track
  LIBS base parameter ${math_libs} container device psi
  SOURCES ${state_track_test_sources}
)

set(cal_edm_test_sources test_cal_edm.cpp ../cal_edm.cpp ../utils/lr_util.cpp)
if(ENABLE_MPI)
  list(APPEND cal_edm_test_sources ../../../source_base/test/mpi_test_main.cpp)
endif()
AddTest(
  TARGET MODULE_LR_cal_edm
  LIBS base parameter ${math_libs} container device psi
  SOURCES ${cal_edm_test_sources}
)
