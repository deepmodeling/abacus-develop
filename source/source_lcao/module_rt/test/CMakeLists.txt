set(SNAP_PSIBETA_CUDA_TEST_SOURCES)
if(USE_CUDA)
  list(APPEND SNAP_PSIBETA_CUDA_TEST_SOURCES
       ../kernels/cuda/snap_psibeta_kernel.cu
       ../kernels/cuda/snap_psibeta_gpu.cu)
endif()

add_library(tddft_test_lib tddft_test.cpp)
target_link_libraries(tddft_test_lib PRIVATE Threads::Threads GTest::gtest_main GTest::gmock_main)
#target_include_directories(tddft_test_lib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${GTEST_INCLUDE_DIRS}>)

AddTest(
  TARGET MODULE_LCAO_tddft_middle_hamilt_test
  LIBS parameter base device tddft_test_lib
  SOURCES middle_hamilt_test.cpp ../middle_hamilt.cpp 
)

AddTest(
  TARGET MODULE_LCAO_tddft_band_energy_test
  LIBS parameter base device tddft_test_lib
  SOURCES band_energy_test.cpp ../band_energy.cpp ../../../source_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET MODULE_LCAO_tddft_norm_psi_test
  LIBS parameter base device tddft_test_lib
  SOURCES norm_psi_test.cpp ../norm_psi.cpp 
)

AddTest(
  TARGET MODULE_LCAO_tddft_upsi_test
  LIBS parameter base device tddft_test_lib
  SOURCES upsi_test1.cpp upsi_test2.cpp upsi_test3.cpp ../upsi.cpp 
)

AddTest(
  TARGET MODULE_LCAO_tddft_propagator_test
  LIBS parameter base device tddft_test_lib
  SOURCES propagator_test1.cpp propagator_test2.cpp propagator_test3.cpp ../propagator.cpp ../propagator_cn2.cpp ../propagator_taylor.cpp ../propagator_etrs.cpp
)

AddTest(
  TARGET MODULE_LCAO_tddft_snap_psibeta_half_test
  LIBS parameter base device orb numerical_atomic_orbitals tddft_test_lib
  SOURCES snap_psibeta_half_tddft_test.cpp ../radial_interpolation.cpp ../snap_proj_half_tddft.cpp ../snap_psb_half_tddft.cpp
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
          ../../../source_io/module_hs/cal_r_overlap_r.cpp
          ../../../source_io/module_hs/single_r_io.cpp
          ../../../source_io/module_hs/rr_sparse_writer.cpp
          ../../../source_pw/module_pwdft/soc.cpp
          ${SNAP_PSIBETA_CUDA_TEST_SOURCES}
)

AddTest(
  TARGET MODULE_LCAO_tddft_radial_interpolation_test
  LIBS base device
  SOURCES radial_interpolation_test.cpp ../radial_interpolation.cpp
)

if(USE_CUDA)
  AddTest(
    TARGET MODULE_LCAO_tddft_radial_interpolation_cuda_test
    LIBS base device
    SOURCES radial_interpolation_cuda_test.cu ../radial_interpolation.cpp
  )
endif()
