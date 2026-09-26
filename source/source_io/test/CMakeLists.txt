abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

AddTest(
  TARGET MODULE_IO_write_bands
  LIBS parameter base device symmetry
  SOURCES write_bands_test.cpp ../module_energy/write_bands.cpp
    ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp
    ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

if(NOT ENABLE_MPI)
  return()
endif()

add_test(NAME MODULE_IO_write_bands_parallel
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:MODULE_IO_write_bands> ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(MODULE_IO_write_bands_parallel PROPERTIES PROCESSORS 4 TIMEOUT 120)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
configure_file(INPUTs ${CMAKE_CURRENT_BINARY_DIR}/INPUTs COPYONLY)

AddTest(
  TARGET MODULE_IO_input_test_para
  LIBS parameter  base device io_input
  SOURCES read_input_ptest.cpp
)

add_test(NAME MODULE_IO_input_test_para_4
      COMMAND mpirun -np 4 ./MODULE_IO_input_test_para
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_read_exit_file_test
  LIBS parameter  base device
  SOURCES read_exit_file_test.cpp
)

add_test(NAME MODULE_IO_read_exit_file_test_para_4
      COMMAND mpirun -np 4 ./MODULE_IO_read_exit_file_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_output_test
  LIBS parameter  base device
  SOURCES output_test.cpp
)

AddTest(
  TARGET MODULE_IO_td_efield_io_test
  LIBS parameter base device
  SOURCES td_efield_io_test.cpp ../module_efield/td_efield_io.cpp
    ../../source_estate/module_pot/td_field.cpp ../../source_estate/module_pot/td_field_profiles.cpp
    ../../source_estate/module_pot/td_field_manager.cpp
)

AddTest(
  TARGET MODULE_IO_td_vector_pot_io_test
  LIBS parameter base device
  SOURCES td_vector_pot_io_test.cpp ../module_efield/td_vector_pot_io.cpp
)

AddTest(
  TARGET MODULE_IO_write_eig_occ_test
  LIBS parameter  base device symmetry
  SOURCES write_eig_occ_test.cpp ../module_energy/write_eig_occ.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp
		../../source_cell/cif_io.cpp ../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_IO_cal_dos
  LIBS parameter  base device
  SOURCES cal_dos_test.cpp ../module_dos/cal_dos.cpp
)

AddTest(
  TARGET MODULE_IO_write_dos_pw
  LIBS parameter  base device symmetry
  SOURCES write_dos_pw_test.cpp ../module_dos/cal_dos.cpp ../module_dos/write_dos_pw.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../module_energy/nscf_fermi_surf.cpp ../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_IO_print_info
  LIBS parameter  base device symmetry cell_info
  SOURCES print_info_test.cpp ../module_output/print_info.cpp ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_IO_write_wfc_nao
  LIBS parameter  base psi device
  SOURCES write_wfc_nao_test.cpp ../module_wf/write_wfc_nao.cpp ../../source_basis/module_ao/parallel_orbitals.cpp
)

configure_file(write_wfc_nao_para.sh ${CMAKE_CURRENT_BINARY_DIR}/write_wfc_nao_para.sh COPYONLY)
find_program(BASH bash)
add_test(NAME MODULE_IO_write_wfc_nao_para
      COMMAND ${BASH} write_wfc_nao_para.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_write_orb_info
  LIBS parameter  base device cell_info
  SOURCES write_orb_info_test.cpp ../module_output/write_orb_info.cpp
)

AddTest(
  TARGET MODULE_IO_parse_args
  LIBS parameter  base device io_input
  SOURCES parse_args_test.cpp ../parse_args.cpp ../input_help.cpp
)

AddTest(
  TARGET MODULE_IO_input_help_test
  LIBS parameter  base device io_input
  SOURCES input_help_test.cpp ../input_help.cpp
)

AddTest(
  TARGET MODULE_IO_bessel_basis_test
  LIBS parameter  base device
  SOURCES bessel_basis_test.cpp ../module_bessel/bessel_basis.cpp
)

AddTest(
  TARGET MODULE_IO_output_log_test
  LIBS parameter base device
  SOURCES ../module_output/output_log.cpp outputlog_test.cpp ../../source_basis/module_pw/test/test_tool.cpp
)

if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_IO_to_qo_test
  LIBS parameter  base device numerical_atomic_orbitals container orb
  SOURCES
  to_qo_test.cpp
  ../module_qo/to_qo_kernel.cpp
  ../module_qo/to_qo_mpi.cpp
  ../module_qo/to_qo_structures.cpp
  ../../source_cell/atom_spec.cpp
  ../../source_cell/parallel_kpoints.cpp
  ../../source_cell/test/support/mock_unitcell.cpp
  ../../source_lcao/center2orb.cpp
)
endif()

AddTest(
  TARGET MODULE_IO_read_wfc_pw_test
  LIBS parameter  base device planewave
  SOURCES read_wfc_pw_test.cpp ../module_wf/read_wfc_pw.cpp ../../source_basis/module_pw/test/test_tool.cpp
)

add_test(NAME MODULE_IO_read_wfc_pw_test_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_read_wfc_pw_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_read_wf2rho_pw_test
  LIBS parameter base device planewave psi symmetry
  SOURCES read_wf2rho_pw_test.cpp ../module_wf/read_wfc_pw.cpp ../module_wf/read_wf2rho_pw.cpp ../../source_basis/module_pw/test/test_tool.cpp  ../../source_estate/module_charge/chg_parallel.cpp ../module_wf/write_wfc_pw.cpp
  ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

add_test(NAME MODULE_IO_read_wf2rho_pw_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_read_wf2rho_pw_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)


if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_IO_numerical_basis_test
  LIBS parameter  base device numerical_atomic_orbitals container orb
  SOURCES numerical_basis_test.cpp
          ../module_bessel/numerical_basis_jyjy.cpp
          ../../source_lcao/center2orb.cpp
)
endif()


AddTest(
  TARGET MODULE_IO_mulliken_test
  LIBS parameter  base device
  SOURCES output_mulliken_test.cpp output_mulliken_mock.cpp ../module_mulliken/output_mulliken.cpp
  ../../source_cell/cell_index.cpp
  ../../source_basis/module_ao/parallel_orbitals.cpp
)

#if(ENABLE_LCAO)
#AddTest(
#  TARGET MODULE_IO_read_wfc_lcao_test
#  LIBS parameter  base device
#  SOURCES read_wfc_lcao_test.cpp ../read_wfc_lcao.cpp
#)

#add_test(NAME MODULE_IO_read_wfc_lcao_test_parallel
#      COMMAND mpirun -np 4 ./MODULE_IO_read_wfc_lcao_test
#      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#)

#endif()

AddTest(
  TARGET MODULE_IO_orb_io_test
  LIBS parameter base device
  SOURCES orb_io_test.cpp
)

add_test(NAME MODULE_IO_orb_io_test_parallel
      COMMAND mpirun -np 4 ./MODULE_IO_orb_io_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_write_dmk
  LIBS parameter base device cell_info symmetry
  SOURCES ../module_dm/test/write_dmk_test.cpp ../module_dm/write_dmk.cpp ../../source_cell/ucell_io.cpp
  ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

add_test(
  NAME MODULE_IO_write_dmk_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_write_dmk
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_IO_read_wfc_nao_test
  LIBS parameter base device
  SOURCES read_wfc_nao_test.cpp ../module_wf/read_wfc_nao.cpp ../module_wf/write_wfc_nao.cpp ../../source_psi/psi.cpp ../../source_basis/module_ao/parallel_orbitals.cpp
)

add_test(
  NAME MODULE_IO_read_wfc_nao_test_parallel
  COMMAND mpirun -np 2 ./MODULE_IO_read_wfc_nao_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

if(ENABLE_LIBRI)
AddTest(
  TARGET MODULE_IO_restart_exx_csr_test
  LIBS parameter base device
  SOURCES
  restart_exx_csr_test.cpp
  tmp_mocks.cpp
)
endif()

AddTest(
  TARGET MODULE_IO_write_elf_logic_test
  SOURCES write_elf_logic_test.cpp
)
