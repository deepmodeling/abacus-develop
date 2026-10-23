set(SNAP_PSIBETA_CUDA_TEST_SOURCES)
if(USE_CUDA)
  list(APPEND SNAP_PSIBETA_CUDA_TEST_SOURCES
       ../kernels/cuda/snap_psibeta_kernel.cu
       ../kernels/cuda/snap_psibeta_gpu.cu)
endif()

add_library(tddft_test_lib tddft_test.cpp)
target_link_libraries(tddft_test_lib PRIVATE Threads::Threads GTest::gtest_main
                                             GTest::gmock_main)
#target_include_directories(tddft_test_lib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${GTEST_INCLUDE_DIRS}>)

# This regression exercises CPU and CPU Tensor-LAPACK propagation. A cuBLASMp
# target also needs its GPU kernels and distributed-library link dependencies.
if(NOT ENABLE_CUBLASMP)
  addtest(
    TARGET
    MODULE_LCAO_tddft_evolve_psi_test
    LIBS
    parameter
    base
    device
    tddft_test_lib
    SOURCES
    test_evolve_psi.cpp
    ../evolve_psi.cpp
    ../middle_hamilt.cpp
    ../norm_psi.cpp
    ../band_energy.cpp
    ../upsi.cpp
    ../solve_propagation.cpp
    ../propagator.cpp
    ../propagator_cn2.cpp
    ../propagator_taylor.cpp
    ../propagator_etrs.cpp
    ../../../source_basis/module_ao/parallel_orbitals.cpp
    ../../../source_io/module_energy/td_eij_io.cpp)
endif()

addtest(
  TARGET
  MODULE_LCAO_tddft_middle_hamilt_test
  LIBS
  parameter
  base
  device
  tddft_test_lib
  SOURCES
  middle_hamilt_test.cpp
  ../middle_hamilt.cpp)

addtest(
  TARGET
  MODULE_LCAO_tddft_band_energy_test
  LIBS
  parameter
  base
  device
  tddft_test_lib
  SOURCES
  band_energy_test.cpp
  ../band_energy.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_io/module_energy/td_eij_io.cpp)

addtest(
  TARGET
  MODULE_LCAO_tddft_norm_psi_test
  LIBS
  parameter
  base
  device
  tddft_test_lib
  SOURCES
  norm_psi_test.cpp
  ../norm_psi.cpp)

addtest(
  TARGET
  MODULE_LCAO_tddft_upsi_test
  LIBS
  parameter
  base
  device
  tddft_test_lib
  SOURCES
  upsi_test1.cpp
  upsi_test2.cpp
  upsi_test3.cpp
  ../upsi.cpp)

addtest(
  TARGET
  MODULE_LCAO_tddft_propagator_test
  LIBS
  parameter
  base
  device
  tddft_test_lib
  SOURCES
  propagator_test1.cpp
  propagator_test2.cpp
  propagator_test3.cpp
  ../propagator.cpp
  ../propagator_cn2.cpp
  ../propagator_taylor.cpp
  ../propagator_etrs.cpp)

addtest(
  TARGET
  MODULE_LCAO_tddft_snap_psibeta_half_test
  LIBS
  parameter
  base
  device
  orb
  numerical_atomic_orbitals
  nao_smoothgen
  tddft_test_lib
  SOURCES
  snap_psb_half_tddft_test.cpp
  ../radial_interpolation.cpp
  ../snap_proj_half_tddft.cpp
  ../snap_psb_half_tddft.cpp
  ../../center2orb.cpp
  ../../center2orb_orb11.cpp
  ../../center2orb_orb21.cpp
  ../../setup_nonlocal.cpp
  ../../../source_cell/atom_spec.cpp
  ../../../source_cell/atom_pseudo.cpp
  ../../../source_cell/pseudo.cpp
  ../../../source_cell/read_pp.cpp
  ../../../source_cell/read_pp_complete.cpp
  ../../../source_cell/read_pp_upf201.cpp
  ../../../source_cell/read_pp_upf100.cpp
  ../../../source_cell/read_pp_vwr.cpp
  ../../../source_cell/read_pp_blps.cpp
  ../../../source_io/module_hs/pos_op_mat.cpp
  ../../../source_io/module_hs/pos_op_basis.cpp
  ../../../source_io/module_hs/pos_op_calc.cpp
  ../../../source_io/module_hs/pos_op_writer.cpp
  ../../../source_io/module_hs/lat_r_csr.cpp
  ../../../source_io/module_hs/pos_op_csr.cpp
  ../../../source_pw/module_pwdft/soc.cpp
  ${SNAP_PSIBETA_CUDA_TEST_SOURCES})

addtest(
  TARGET
  MODULE_LCAO_tddft_rad_interp_test
  LIBS
  base
  device
  SOURCES
  rad_interp_test.cpp
  ../radial_interpolation.cpp)

if(USE_CUDA)
  addtest(
    TARGET
    MODULE_LCAO_tddft_rad_interp_cuda_test
    LIBS
    base
    device
    SOURCES
    rad_interp_cuda_test.cu
    ../radial_interpolation.cpp)
endif()

addtest(
  TARGET MODULE_LCAO_tddft_boundary_fix_test
  LIBS parameter base device container cell symmetry neighbor neighbor_search
       domain_decomposition psi
  SOURCES
  test_boundary_fix.cpp
  ../boundary_fix.cpp
  ../../../source_base/test/mpi_test_main.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_md/md_base.cpp
  ../../../source_md/md_func.cpp
  ../../../source_io/module_output/output_log.cpp
  ../../../source_io/module_output/print_info.cpp)
add_test(
  NAME MODULE_LCAO_tddft_boundary_fix_test_np4
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
          $<TARGET_FILE:MODULE_LCAO_tddft_boundary_fix_test>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(
  MODULE_LCAO_tddft_boundary_fix_test MODULE_LCAO_tddft_boundary_fix_test_np4
  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;I_MPI_PIN=0" TIMEOUT 60)
