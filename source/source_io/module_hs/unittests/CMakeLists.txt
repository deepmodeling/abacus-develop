AddTest(
  TARGET MODULE_IO_hsk_writer_test
  LIBS parameter base device
  SOURCES
  test_hsk_writer.cpp
  ../hsk_writer.cpp
  ../hs_dense_io.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
)

add_test(
  NAME MODULE_IO_hsk_writer_test_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_hsk_writer_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_lat_r_csr_test
  LIBS parameter
  SOURCES test_lat_r_csr.cpp ../lat_r_csr.cpp
	../../../source_base/global_variable.cpp
	../../../source_base/parallel_reduce.cpp
	../../../source_base/parallel_common.cpp
	../../../source_base/parallel_global.cpp
  ../../../source_base/parallel_comm.cpp
  ../../../source_base/tool_quit.cpp ../../../source_base/global_file.cpp ../../../source_base/global_function.cpp ../../../source_base/memory_recorder.cpp ../../../source_base/timer.cpp
)

if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_IO_vxc_op_tools_test
  LIBS parameter base device
  SOURCES
  test_vxc_op_tools.cpp
  ../vxc_op_tools.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET MODULE_IO_hs_dense_io_test
  LIBS parameter
  SOURCES
  test_hs_dense_io.cpp
  ../hs_dense_io.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_base/global_variable.cpp
  ../../../source_base/parallel_2d.cpp
  ../../../source_base/parallel_reduce.cpp
  ../../../source_base/parallel_common.cpp
  ../../../source_base/parallel_global.cpp
  ../../../source_base/parallel_comm.cpp
  ../../../source_base/tool_quit.cpp
  ../../../source_base/tool_title.cpp
  ../../../source_base/global_file.cpp
  ../../../source_base/global_function.cpp
  ../../../source_base/memory_recorder.cpp
  ../../../source_base/timer.cpp
)

add_test(
  NAME MODULE_IO_hs_dense_io_test_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_hs_dense_io_test
          --gtest_filter=HsDenseIoSaveMat.BinaryGammaUpperTriangle:HsDenseIoSaveMat.BinaryComplexUpperTriangleMpiReduction:HsDenseIoSaveMat.BinaryAppendCompleteRecordsAndOverwrite
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_hsr_writer_test
  LIBS parameter base device hcontainer
  SOURCES
  test_hsr_writer.cpp
  ../hsr_writer.cpp
  ../../../source_cell/ucell_io.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_hamilt/module_hcontainer/test/tmp_mocks.cpp
)

add_test(
  NAME MODULE_IO_hsr_writer_test_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_hsr_writer_test
          --gtest_filter=HsrWriterIo.HContainerBinaryMpiGatherWritesCompleteFiles
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_angmom_op_mat_test
  LIBS parameter base device neighbor
  SOURCES
  test_angmom_op_mat.cpp
  ../angmom_op_mat.cpp
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
  TARGET MODULE_IO_pos_op_csr_test
  LIBS parameter base device
  SOURCES
  test_pos_op_csr.cpp
  ../pos_op_csr.cpp
)

AddTest(
  TARGET MODULE_IO_hs_sparse_io_test
  LIBS parameter base device hcontainer
  SOURCES
  test_hs_sparse_io.cpp
  ../hs_sparse_io.cpp
  ../lat_r_csr.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_hamilt/module_hcontainer/test/tmp_mocks.cpp
)

AddTest(
  TARGET MODULE_IO_dhs_sparse_writer_test
  LIBS parameter base device hcontainer
  SOURCES
  test_dhs_sparse_writer.cpp
  ../dhs_sparse_writer.cpp
  ../hs_sparse_io.cpp
  ../lat_r_csr.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_hamilt/module_hcontainer/test/tmp_mocks.cpp
)

AddTest(
  TARGET MODULE_IO_output_mat_sparse_test
  LIBS parameter base device
  SOURCES
  test_output_mat_sparse.cpp
)
endif()
