set(SNAP_PSIBETA_CUDA_TEST_SOURCES)
if(USE_CUDA)
  list(APPEND SNAP_PSIBETA_CUDA_TEST_SOURCES
       ../kernels/cuda/snap_psibeta_kernel.cu
       ../kernels/cuda/snap_psibeta_gpu.cu)
endif()

addtest(
  TARGET MODULE_LCAO_tddft_snap_psibeta_half_test
  LIBS parameter
       base
       device
       orb
       numerical_atomic_orbitals
       nao_smoothgen
       tddft_test_lib
  SOURCES test_snap_psb_half_tddft.cpp
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
          ../../../source_cell/ucell_io.cpp
          ../../../source_lcao/module_op/pos_op_mat.cpp
          ../../../source_lcao/module_op/pos_op_basis.cpp
          ../../../source_lcao/module_op/pos_op_calc.cpp
          ../../../source_lcao/module_op/pos_op_writer.cpp
          ../../../source_io/module_hs/lat_r_csr.cpp
          ../../../source_lcao/module_op/pos_op_csr.cpp
          ../../../source_cell/soc.cpp
          ${SNAP_PSIBETA_CUDA_TEST_SOURCES})
