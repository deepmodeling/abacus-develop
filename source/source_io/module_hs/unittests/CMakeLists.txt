AddTest(
  TARGET MODULE_IO_write_hsk_binary_test
  LIBS parameter base device
  SOURCES write_hsk_binary_test.cpp ../../../source_basis/module_ao/parallel_orbitals.cpp
)

add_test(
  NAME MODULE_IO_write_hsk_binary_test_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_write_hsk_binary_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_single_R_test
  LIBS parameter
  SOURCES single_r_io_test.cpp ../single_r_io.cpp
	../../../source_base/global_variable.cpp
	../../../source_base/parallel_reduce.cpp
	../../../source_base/parallel_common.cpp
	../../../source_base/parallel_global.cpp
  ../../../source_base/parallel_comm.cpp
  ../../../source_base/tool_quit.cpp ../../../source_base/global_file.cpp ../../../source_base/global_function.cpp ../../../source_base/memory_recorder.cpp ../../../source_base/timer.cpp
)

if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_IO_cal_plpr_test
  LIBS parameter base device neighbor
  SOURCES
  cal_plpr_test.cpp
  ../cal_plpr.cpp
  ../../../source_basis/module_ao/orb_atomic_lm.cpp
  ../../../source_basis/module_ao/orb_atomic.cpp
  ../../../source_basis/module_nao/radial_set.cpp
  ../../../source_basis/module_nao/numerical_radial.cpp
  ../../../source_basis/module_nao/beta_radials.cpp
  ../../../source_basis/module_nao/hydrogen_radials.cpp
  ../../../source_basis/module_nao/pswfc_radials.cpp
  ../../../source_basis/module_nao/atomic_radials.cpp
  ../../../source_basis/module_nao/sphbes_radials.cpp
  ../../../source_basis/module_nao/two_center_integrator.cpp
  ../../../source_basis/module_nao/two_center_table.cpp
  ../../../source_basis/module_nao/real_gaunt_table.cpp
  ../../../source_basis/module_nao/radial_collection.cpp
)

AddTest(
  TARGET MODULE_IO_write_hs_r_compat_test
  LIBS parameter base device hcontainer
  SOURCES
  write_hs_r_compat_test.cpp
  ../write_hs_r.cpp
  ../write_hs_sparse.cpp
  ../single_r_io.cpp
  ../rr_sparse_writer.cpp
  ../../module_dm/write_dmr.cpp
  ../../../source_cell/ucell_io.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_hamilt/module_hcontainer/test/tmp_mocks.cpp
)

add_test(
  NAME MODULE_IO_write_hsr_binary_test_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_write_hs_r_compat_test
          --gtest_filter=WriteHsRCompatibility.HContainerBinaryMpiGatherWritesCompleteFiles
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
endif()
