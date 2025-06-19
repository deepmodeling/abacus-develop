remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES INPUTs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_IO_input_test_para
  LIBS parameter  ${math_libs} base device io_input 
  SOURCES read_input_ptest.cpp
)

add_test(NAME MODULE_IO_input_test_para_4
      COMMAND mpirun -np 4 ./MODULE_IO_input_test_para
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_read_exit_file_test
  LIBS parameter  ${math_libs} base device 
  SOURCES read_exit_file_test.cpp ../read_exit_file.cpp
)

add_test(NAME MODULE_IO_read_exit_file_test_para_4
      COMMAND mpirun -np 4 ./MODULE_IO_read_exit_file_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_winput_test
  LIBS parameter ${math_libs} base device 
  SOURCES winput_test.cpp ../winput.cpp
)

AddTest(
  TARGET MODULE_IO_output_test
  LIBS parameter  ${math_libs} base device
  SOURCES output_test.cpp ../output.cpp
)

AddTest(
  TARGET MODULE_IO_binstream_test
  SOURCES binstream_test.cpp ../binstream.cpp
)

AddTest(
  TARGET MODULE_IO_write_istate_info_test
  LIBS parameter  ${math_libs} base device symmetry 
  SOURCES write_istate_info_test.cpp ../write_istate_info.cpp ../output.cpp ../../module_cell/parallel_kpoints.cpp ../../module_cell/klist.cpp ../../module_cell/k_vector_utils.cpp
		../cif_io.cpp
)

AddTest(
  TARGET MODULE_IO_cal_dos
  LIBS parameter  ${math_libs} base device
  SOURCES cal_dos_test.cpp ../cal_dos.cpp
)

AddTest(
  TARGET MODULE_IO_write_dos_pw
  LIBS parameter  ${math_libs} base device symmetry
  SOURCES write_dos_pw_test.cpp ../cal_dos.cpp ../write_dos_pw.cpp ../output.cpp ../../module_cell/parallel_kpoints.cpp ../../module_cell/klist.cpp ../nscf_fermi_surf.cpp ../../module_cell/k_vector_utils.cpp
)

AddTest(
  TARGET MODULE_IO_print_info
  LIBS parameter  ${math_libs} base device symmetry cell_info 
  SOURCES print_info_test.cpp ../print_info.cpp ../output.cpp ../../module_cell/klist.cpp ../../module_cell/parallel_kpoints.cpp ../../module_cell/k_vector_utils.cpp
)

AddTest(
  TARGET MODULE_IO_single_R_test
  LIBS parameter  ${math_libs}
  SOURCES single_R_io_test.cpp ../single_R_io.cpp
  	../../module_base/global_variable.cpp
	../../module_base/parallel_reduce.cpp
	../../module_base/parallel_common.cpp
	../../module_base/parallel_global.cpp
  ../../module_base/parallel_comm.cpp
)

AddTest(
  TARGET MODULE_IO_write_wfc_nao
  LIBS parameter  ${math_libs} base psi device
  SOURCES write_wfc_nao_test.cpp ../filename.cpp ../write_wfc_nao.cpp ../../module_basis/module_ao/parallel_orbitals.cpp ../binstream.cpp
)

install(FILES write_wfc_nao_para.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME MODULE_IO_write_wfc_nao_para
      COMMAND ${BASH} write_wfc_nao_para.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_write_orb_info
  LIBS parameter  ${math_libs} base device cell_info 
  SOURCES write_orb_info_test.cpp ../write_orb_info.cpp ../output.cpp
)

AddTest(
  TARGET MODULE_IO_parse_args
  SOURCES parse_args_test.cpp ../parse_args.cpp
)

AddTest(
  TARGET MODULE_IO_bessel_basis_test
  LIBS parameter  ${math_libs} base device
  SOURCES bessel_basis_test.cpp ../bessel_basis.cpp
)

AddTest(
  TARGET MODULE_IO_output_log_test
  LIBS parameter base ${math_libs} device 
  SOURCES ../output_log.cpp outputlog_test.cpp ../../module_basis/module_pw/test/test_tool.cpp
)

AddTest(
  TARGET MODULE_IO_sparse_matrix_test
  LIBS parameter  base ${math_libs} device
  SOURCES sparse_matrix_test.cpp ../sparse_matrix.cpp
)

AddTest(
  TARGET MODULE_IO_file_reader_test
  LIBS parameter  base ${math_libs} device
  SOURCES file_reader_test.cpp ../file_reader.cpp
)

AddTest(
  TARGET MODULE_IO_csr_reader_test
  LIBS parameter  base ${math_libs} device
  SOURCES csr_reader_test.cpp ../csr_reader.cpp ../file_reader.cpp ../sparse_matrix.cpp
)

AddTest(
  TARGET MODULE_IO_read_rhog_test
  LIBS parameter  ${math_libs} base device planewave
  SOURCES read_rhog_test.cpp ../rhog_io.cpp ../binstream.cpp ../../module_basis/module_pw/test/test_tool.cpp
)

if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_IO_to_qo_test
  LIBS parameter  base ${math_libs} device numerical_atomic_orbitals container orb 
  SOURCES 
  to_qo_test.cpp
  ../to_qo_kernel.cpp
  ../to_qo_mpi.cpp
  ../to_qo_structures.cpp
  ../../module_cell/atom_spec.cpp
  ../../module_cell/parallel_kpoints.cpp
  ../../module_cell/test/support/mock_unitcell.cpp
  ../../module_hamilt_lcao/hamilt_lcaodft/center2_orb.cpp
  ../orb_io.cpp
)
endif()

AddTest(
  TARGET MODULE_IO_read_wfc_pw_test
  LIBS parameter  base ${math_libs} device planewave
  SOURCES read_wfc_pw_test.cpp ../read_wfc_pw.cpp ../binstream.cpp ../../module_basis/module_pw/test/test_tool.cpp
)

add_test(NAME MODULE_IO_read_wfc_pw_test_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_read_wfc_pw_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_read_wf2rho_pw_test
  LIBS parameter base ${math_libs} device planewave psi
  SOURCES read_wf2rho_pw_test.cpp ../read_wfc_pw.cpp ../read_wf2rho_pw.cpp ../binstream.cpp ../../module_basis/module_pw/test/test_tool.cpp  ../../module_elecstate/module_charge/charge_mpi.cpp ../filename.cpp ../write_wfc_pw.cpp
)

add_test(NAME MODULE_IO_read_wf2rho_pw_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_read_wf2rho_pw_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)


AddTest(
  TARGET MODULE_IO_numerical_basis_test
  LIBS parameter  base ${math_libs} device numerical_atomic_orbitals container orb 
  SOURCES numerical_basis_test.cpp 
          ../numerical_basis_jyjy.cpp 
          ../../module_hamilt_lcao/hamilt_lcaodft/center2_orb.cpp
          ../orb_io.cpp
)


AddTest(
  TARGET MODULE_IO_mulliken_test
  LIBS parameter  base ${math_libs} device
  SOURCES output_mulliken_test.cpp output_mulliken_mock.cpp ../output_mulliken.cpp
  ../../module_cell/cell_index.cpp
  ../../module_basis/module_ao/parallel_orbitals.cpp
  ../orb_io.cpp
)

#if(ENABLE_LCAO)
#AddTest(
#  TARGET MODULE_IO_read_wfc_lcao_test
#  LIBS parameter  base ${math_libs} device
#  SOURCES read_wfc_lcao_test.cpp ../read_wfc_lcao.cpp
#)

#add_test(NAME MODULE_IO_read_wfc_lcao_test_parallel
#      COMMAND mpirun -np 4 ./MODULE_IO_read_wfc_lcao_test
#      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#)

#endif()

AddTest(
  TARGET MODULE_IO_cif_io_test
  LIBS parameter base ${math_libs} device
  SOURCES cif_io_test.cpp ../cif_io.cpp
)

add_test(NAME MODULE_IO_cif_io_test_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_cif_io_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_orb_io_test
  LIBS parameter base ${math_libs} device
  SOURCES orb_io_test.cpp ../orb_io.cpp
)

add_test(NAME MODULE_IO_orb_io_test_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_orb_io_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_dmk_io
  LIBS parameter ${math_libs} base device cell_info 
  SOURCES io_dmk_test.cpp ../io_dmk.cpp ../output.cpp
)

add_test(
  NAME MODULE_IO_dmk_io_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_dmk_io
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_read_wfc_nao_test
  LIBS parameter ${math_libs} base device
  SOURCES read_wfc_nao_test.cpp ../read_wfc_nao.cpp ../../module_psi/psi.cpp ../../module_basis/module_ao/parallel_orbitals.cpp
)

add_test(
  NAME MODULE_IO_read_wfc_nao_test_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_read_wfc_nao_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_IO_cal_pLpR_test
  LIBS parameter base ${math_libs} device neighbor
  SOURCES 
  cal_pLpR_test.cpp 
  ../cal_pLpR.cpp
  ../../module_basis/module_ao/ORB_atomic_lm.cpp
  ../../module_basis/module_ao/ORB_atomic.cpp
  ../../module_basis/module_nao/radial_set.cpp
  ../../module_basis/module_nao/numerical_radial.cpp
  ../../module_basis/module_nao/beta_radials.cpp
  ../../module_basis/module_nao/hydrogen_radials.cpp
  ../../module_basis/module_nao/pswfc_radials.cpp
  ../../module_basis/module_nao/atomic_radials.cpp
  ../../module_basis/module_nao/sphbes_radials.cpp
  ../../module_basis/module_nao/two_center_integrator.cpp
  ../../module_basis/module_nao/two_center_table.cpp
  ../../module_basis/module_nao/real_gaunt_table.cpp
  ../../module_basis/module_nao/radial_collection.cpp
)
endif()
