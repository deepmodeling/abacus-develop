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
  SOURCES read_exit_file_test.cpp ../module_output/read_exit_file.cpp
)

add_test(NAME MODULE_IO_read_exit_file_test_para_4
      COMMAND mpirun -np 4 ./MODULE_IO_read_exit_file_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_output_test
  LIBS parameter  ${math_libs} base device
  SOURCES output_test.cpp ../module_output/output.cpp
)

AddTest(
  TARGET MODULE_IO_binstream_test
  SOURCES binstream_test.cpp ../module_output/binstream.cpp
)

AddTest(
  TARGET MODULE_IO_write_eig_occ_test
  LIBS parameter  ${math_libs} base device symmetry
  SOURCES write_eig_occ_test.cpp ../module_energy/write_eig_occ.cpp ../module_output/output.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/klist.cpp ../../source_cell/k_vector_utils.cpp
		../module_output/cif_io.cpp
)

AddTest(
  TARGET MODULE_IO_cal_dos
  LIBS parameter  ${math_libs} base device
  SOURCES cal_dos_test.cpp ../module_dos/cal_dos.cpp
)

AddTest(
  TARGET MODULE_IO_write_dos_pw
  LIBS parameter  ${math_libs} base device symmetry
  SOURCES write_dos_pw_test.cpp ../module_dos/cal_dos.cpp ../module_dos/write_dos_pw.cpp ../module_output/output.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/klist.cpp ../module_energy/nscf_fermi_surf.cpp ../../source_cell/k_vector_utils.cpp
)

AddTest(
  TARGET MODULE_IO_print_info
  LIBS parameter  ${math_libs} base device symmetry cell_info
  SOURCES print_info_test.cpp ../module_output/print_info.cpp ../module_output/output.cpp ../../source_cell/klist.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/k_vector_utils.cpp
)

AddTest(
  TARGET MODULE_IO_single_R_test
  LIBS parameter  ${math_libs}
  SOURCES single_R_io_test.cpp ../module_hs/single_R_io.cpp
	../../source_base/global_variable.cpp
	../../source_base/parallel_reduce.cpp
	../../source_base/parallel_common.cpp
	../../source_base/parallel_global.cpp
  ../../source_base/parallel_comm.cpp
  ../../source_base/tool_quit.cpp ../../source_base/global_file.cpp ../../source_base/global_function.cpp ../../source_base/memory.cpp ../../source_base/timer.cpp
)

AddTest(
  TARGET MODULE_IO_write_wfc_nao
  LIBS parameter  ${math_libs} base psi device
  SOURCES write_wfc_nao_test.cpp ../module_output/filename.cpp ../module_wf/write_wfc_nao.cpp ../../source_basis/module_ao/parallel_orbitals.cpp ../module_output/binstream.cpp
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
  SOURCES write_orb_info_test.cpp ../module_output/write_orb_info.cpp ../module_output/output.cpp
)

AddTest(
  TARGET MODULE_IO_parse_args
  LIBS parameter  ${math_libs} base device io_input
  SOURCES parse_args_test.cpp ../parse_args.cpp
)

AddTest(
  TARGET MODULE_IO_bessel_basis_test
  LIBS parameter  ${math_libs} base device
  SOURCES bessel_basis_test.cpp ../module_bessel/bessel_basis.cpp
)

AddTest(
  TARGET MODULE_IO_output_log_test
  LIBS parameter base ${math_libs} device
  SOURCES ../module_output/output_log.cpp outputlog_test.cpp ../../source_basis/module_pw/test/test_tool.cpp
)

AddTest(
  TARGET MODULE_IO_sparse_matrix_test
  LIBS parameter  base ${math_libs} device
  SOURCES sparse_matrix_test.cpp ../module_output/sparse_matrix.cpp
)

AddTest(
  TARGET MODULE_IO_file_reader_test
  LIBS parameter  base ${math_libs} device
  SOURCES file_reader_test.cpp ../module_output/file_reader.cpp
)

AddTest(
  TARGET MODULE_IO_csr_reader_test
  LIBS parameter  base ${math_libs} device
  SOURCES csr_reader_test.cpp ../module_output/csr_reader.cpp ../module_output/file_reader.cpp ../module_output/sparse_matrix.cpp
)

AddTest(
  TARGET MODULE_IO_read_rhog_test
  LIBS parameter  ${math_libs} base device planewave
  SOURCES read_rhog_test.cpp ../module_chgpot/rhog_io.cpp ../module_output/binstream.cpp ../../source_basis/module_pw/test/test_tool.cpp
)

if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_IO_to_qo_test
  LIBS parameter  base ${math_libs} device numerical_atomic_orbitals container orb
  SOURCES
  to_qo_test.cpp
  ../module_qo/to_qo_kernel.cpp
  ../module_qo/to_qo_mpi.cpp
  ../module_qo/to_qo_structures.cpp
  ../../source_cell/atom_spec.cpp
  ../../source_cell/parallel_kpoints.cpp
  ../../source_cell/test/support/mock_unitcell.cpp
  ../../source_lcao/center2_orb.cpp
  ../module_output/orb_io.cpp
)
endif()

AddTest(
  TARGET MODULE_IO_read_wfc_pw_test
  LIBS parameter  base ${math_libs} device planewave
  SOURCES read_wfc_pw_test.cpp ../module_wf/read_wfc_pw.cpp ../module_output/binstream.cpp ../../source_basis/module_pw/test/test_tool.cpp
)

add_test(NAME MODULE_IO_read_wfc_pw_test_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_read_wfc_pw_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_read_wf2rho_pw_test
  LIBS parameter base ${math_libs} device planewave psi
  SOURCES read_wf2rho_pw_test.cpp ../module_wf/read_wfc_pw.cpp ../module_wf/read_wf2rho_pw.cpp ../module_output/binstream.cpp ../../source_basis/module_pw/test/test_tool.cpp  ../../source_estate/module_charge/charge_mpi.cpp ../module_output/filename.cpp ../module_wf/write_wfc_pw.cpp
)

add_test(NAME MODULE_IO_read_wf2rho_pw_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_read_wf2rho_pw_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)


AddTest(
  TARGET MODULE_IO_numerical_basis_test
  LIBS parameter  base ${math_libs} device numerical_atomic_orbitals container orb
  SOURCES numerical_basis_test.cpp
          ../module_bessel/numerical_basis_jyjy.cpp
          ../../source_lcao/center2_orb.cpp
          ../module_output/orb_io.cpp
)


AddTest(
  TARGET MODULE_IO_mulliken_test
  LIBS parameter  base ${math_libs} device
  SOURCES output_mulliken_test.cpp output_mulliken_mock.cpp ../module_mulliken/output_mulliken.cpp
  ../../source_cell/cell_index.cpp
  ../../source_basis/module_ao/parallel_orbitals.cpp
  ../module_output/orb_io.cpp
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
  SOURCES cif_io_test.cpp ../module_output/cif_io.cpp
)

add_test(NAME MODULE_IO_cif_io_test_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_cif_io_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_orb_io_test
  LIBS parameter base ${math_libs} device
  SOURCES orb_io_test.cpp ../module_output/orb_io.cpp
)

add_test(NAME MODULE_IO_orb_io_test_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_orb_io_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_write_dmk
  LIBS parameter ${math_libs} base device cell_info
  SOURCES ../module_dm/test/write_dmk_test.cpp ../module_dm/write_dmk.cpp ../module_output/output.cpp
)

add_test(
  NAME MODULE_IO_write_dmk_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_write_dmk
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_read_wfc_nao_test
  LIBS parameter ${math_libs} base device
  SOURCES read_wfc_nao_test.cpp ../module_wf/read_wfc_nao.cpp ../../source_psi/psi.cpp ../../source_basis/module_ao/parallel_orbitals.cpp
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
  ../module_hs/cal_pLpR.cpp
  ../../source_basis/module_ao/ORB_atomic_lm.cpp
  ../../source_basis/module_ao/ORB_atomic.cpp
  ../../source_basis/module_nao/radial_set.cpp
  ../../source_basis/module_nao/numerical_radial.cpp
  ../../source_basis/module_nao/beta_radials.cpp
  ../../source_basis/module_nao/hydrogen_radials.cpp
  ../../source_basis/module_nao/pswfc_radials.cpp
  ../../source_basis/module_nao/atomic_radials.cpp
  ../../source_basis/module_nao/sphbes_radials.cpp
  ../../source_basis/module_nao/two_center_integrator.cpp
  ../../source_basis/module_nao/two_center_table.cpp
  ../../source_basis/module_nao/real_gaunt_table.cpp
  ../../source_basis/module_nao/radial_collection.cpp
)
endif()

AddTest(
  TARGET MODULE_IO_write_elf_logic_test
  SOURCES write_elf_logic_test.cpp
)

