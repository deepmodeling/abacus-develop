abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

addtest(
  TARGET MODULE_IO_write_bands
  LIBS parameter base device symmetry
  SOURCES test_write_bands.cpp
          ../module_energy/write_bands.cpp
          ../../source_cell/klist.cpp
          ../../source_cell/klist_io.cpp
          ../../source_cell/parallel_kpoints.cpp
          ../../source_cell/reciprocal_grid.cpp)

addtest(TARGET MODULE_IO_spin_tag SOURCES test_spin_tag.cpp)

if(NOT ENABLE_MPI)
  return()
endif()

add_executable(MODULE_IO_input_error_mpi test_read_input.cpp)
target_link_libraries(
  MODULE_IO_input_error_mpi PRIVATE parameter base device io_input
                                    abacus::linalg_libs)
if(ENABLE_COVERAGE)
  add_coverage(MODULE_IO_input_error_mpi)
endif()
foreach(
  input_case IN
  ITEMS missing_file
        missing_header
        unknown_parameter
        duplicate_parameter
        duplicate_string_eof
        invalid_value
        invalid_dependency
        valid)
  add_test(
    NAME MODULE_IO_input_error_mpi_${input_case}
    COMMAND
      ${CMAKE_COMMAND} "-DREADER=$<TARGET_FILE:MODULE_IO_input_error_mpi>"
      "-DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}"
      "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
      "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
      "-DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}" "-DINPUT_CASE=${input_case}"
      "-DCASE_DIR=${CMAKE_CURRENT_BINARY_DIR}/input_error_${input_case}" -P
      ${CMAKE_CURRENT_SOURCE_DIR}/test_read_input.cmake)
  set_tests_properties(
    MODULE_IO_input_error_mpi_${input_case}
    PROPERTIES PROCESSORS 4 TIMEOUT 45 ENVIRONMENT "OMP_NUM_THREADS=1")
endforeach()

add_test(
  NAME MODULE_IO_write_bands_parallel
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:MODULE_IO_write_bands> ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(MODULE_IO_write_bands_parallel PROPERTIES PROCESSORS 4
                                                               TIMEOUT 120)

abacus_copy_test_files(support INPUTs)

addtest(
  TARGET MODULE_IO_input_test_para
  LIBS parameter base device io_input
  SOURCES test_read_input_parallel.cpp)

add_test(
  NAME MODULE_IO_input_test_para_4
  COMMAND mpirun -np 4 "$<TARGET_FILE:MODULE_IO_input_test_para>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

addtest(
  TARGET MODULE_IO_read_exit_file_test
  LIBS parameter base device
  SOURCES test_read_exit_file.cpp)

add_test(
  NAME MODULE_IO_read_exit_file_test_para_4
  COMMAND mpirun -np 4 "$<TARGET_FILE:MODULE_IO_read_exit_file_test>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

addtest(
  TARGET MODULE_IO_output_test
  LIBS parameter base device
  SOURCES test_output.cpp)

addtest(
  TARGET MODULE_IO_td_efield_io_test
  LIBS parameter base device
  SOURCES test_td_efield_io.cpp
          ../module_efield/td_efield_io.cpp
          ../../source_estate/module_pot/td_field.cpp
          ../../source_estate/module_pot/td_field_profiles.cpp
          ../../source_estate/module_pot/td_field_manager.cpp)

addtest(
  TARGET MODULE_IO_td_vecpot_io_test
  LIBS parameter base device
  SOURCES test_td_vecpot_io.cpp ../module_efield/td_vecpot_io.cpp)

addtest(
  TARGET MODULE_IO_write_eig_occ_test
  LIBS parameter base device symmetry
  SOURCES test_write_eig_occ.cpp
          ../module_energy/write_eig_occ.cpp
          ../../source_cell/parallel_kpoints.cpp
          ../../source_cell/klist.cpp
          ../../source_cell/klist_io.cpp
          ../../source_cell/cif_io.cpp
          ../../source_cell/reciprocal_grid.cpp)

# DOS tests moved to ../module_dos/unittests/

addtest(
  TARGET MODULE_IO_print_info
  LIBS parameter base device symmetry cell_info
  SOURCES test_print_info.cpp
          ../module_output/print_info.cpp
          ../../source_cell/klist.cpp
          ../../source_cell/klist_io.cpp
          ../../source_cell/parallel_kpoints.cpp
          ../../source_cell/reciprocal_grid.cpp)

addtest(
  TARGET MODULE_IO_write_wfc_nao
  LIBS parameter base psi device
  SOURCES test_write_wfc_nao.cpp ../module_wf/write_wfc_nao.cpp
          ../../source_basis/module_ao/parallel_orbitals.cpp)

abacus_copy_test_files(write_wfc_nao_para.sh)
find_program(BASH bash)
add_test(
  NAME MODULE_IO_write_wfc_nao_para
  COMMAND ${BASH} write_wfc_nao_para.sh "$<TARGET_FILE:MODULE_IO_write_wfc_nao>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

addtest(
  TARGET MODULE_IO_write_orb_info
  LIBS parameter base device cell_info
  SOURCES test_write_orb_info.cpp ../module_output/write_orb_info.cpp)

addtest(
  TARGET MODULE_IO_parse_args
  LIBS parameter base device io_input
  SOURCES test_parse_args.cpp ../parse_args.cpp ../input_help.cpp)

addtest(
  TARGET MODULE_IO_test_read_inp_sys
  LIBS base device
  SOURCES test_read_inp_sys.cpp)

addtest(
  TARGET MODULE_IO_input_help_test
  LIBS parameter base device io_input
  SOURCES test_input_help.cpp ../input_help.cpp)

addtest(
  TARGET MODULE_IO_bessel_basis_test
  LIBS parameter base device
  SOURCES test_bessel_basis.cpp ../module_bessel/bessel_basis.cpp)

addtest(
  TARGET MODULE_IO_output_log_test
  LIBS parameter base device
  SOURCES ../module_output/output_log.cpp test_output_log.cpp
          ../../source_basis/module_pw/test/test_tool.cpp)

if(ENABLE_LCAO)
  addtest(
    TARGET MODULE_IO_to_qo_test
    LIBS parameter
         base
         device
         numerical_atomic_orbitals
         nao_smoothgen
         container
         orb
    SOURCES test_to_qo.cpp
            ../module_qo/to_qo_kernel.cpp
            ../module_qo/to_qo_mpi.cpp
            ../module_qo/to_qo_structures.cpp
            ../../source_cell/atom_spec.cpp
            ../../source_cell/parallel_kpoints.cpp
            ../../source_cell/test/support/mock_unitcell.cpp
            ../../source_lcao/center2orb.cpp)
endif()

addtest(
  TARGET MODULE_IO_read_wfc_pw_test
  LIBS parameter base device planewave
  SOURCES test_read_wfc_pw.cpp ../module_wf/read_wfc_pw.cpp
          ../../source_basis/module_pw/test/test_tool.cpp)

add_test(
  NAME MODULE_IO_read_wfc_pw_test_parallel
  COMMAND mpirun -np 4 "$<TARGET_FILE:MODULE_IO_read_wfc_pw_test>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

addtest(
  TARGET MODULE_IO_read_wf2rho_pw_test
  LIBS parameter base device planewave psi symmetry
  SOURCES test_read_wf2rho_pw.cpp
          ../module_wf/read_wfc_pw.cpp
          ../module_wf/read_wf2rho_pw.cpp
          ../../source_basis/module_pw/test/test_tool.cpp
          ../../source_estate/module_charge/charge.cpp
          ../../source_estate/module_charge/chg_tools.cpp
          ../../source_estate/module_charge/chg_parallel.cpp
          ../module_wf/write_wfc_pw.cpp
          ../../source_cell/klist.cpp
          ../../source_cell/klist_io.cpp
          ../../source_cell/parallel_kpoints.cpp
          ../../source_cell/reciprocal_grid.cpp)

add_test(
  NAME MODULE_IO_read_wf2rho_pw_parallel
  COMMAND mpirun -np 4 "$<TARGET_FILE:MODULE_IO_read_wf2rho_pw_test>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_LCAO)
  addtest(
    TARGET MODULE_IO_numerical_basis_test
    LIBS parameter
         base
         device
         numerical_atomic_orbitals
         nao_smoothgen
         container
         orb
    SOURCES test_numerical_basis.cpp ../module_bessel/numerical_basis_jyjy.cpp
            ../../source_lcao/center2orb.cpp)
endif()

addtest(
  TARGET MODULE_IO_mulliken_test
  LIBS parameter base device
  SOURCES test_output_mulliken.cpp
          output_mulliken_mock.cpp
          ../module_mulliken/output_mulliken.cpp
          ../../source_cell/cell_index.cpp
          ../../source_basis/module_ao/parallel_orbitals.cpp)

#if(ENABLE_LCAO)
#AddTest(
#  TARGET MODULE_IO_read_wfc_lcao_test
#  LIBS parameter  base device
#  SOURCES test_read_wfc_lcao.cpp ../read_wfc_lcao.cpp
#)

#add_test(NAME MODULE_IO_read_wfc_lcao_test_parallel
#      COMMAND mpirun -np 4 ./MODULE_IO_read_wfc_lcao_test
#      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#)

#endif()

addtest(
  TARGET MODULE_IO_orb_io_test
  LIBS parameter base device
  SOURCES test_orb_io.cpp)

add_test(
  NAME MODULE_IO_orb_io_test_parallel
  COMMAND mpirun -np 4 "$<TARGET_FILE:MODULE_IO_orb_io_test>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

addtest(
  TARGET MODULE_IO_write_dmk
  LIBS parameter base device cell_info symmetry
  SOURCES ../module_dm/test/write_dmk_test.cpp
          ../module_dm/write_dmk.cpp
          ../../source_cell/ucell_io.cpp
          ../../source_cell/klist.cpp
          ../../source_cell/klist_io.cpp
          ../../source_cell/parallel_kpoints.cpp
          ../../source_cell/reciprocal_grid.cpp)

add_test(
  NAME MODULE_IO_write_dmk_parallel
  COMMAND mpirun -np 2 "$<TARGET_FILE:MODULE_IO_write_dmk>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_LCAO)
  addtest(
    TARGET MODULE_IO_write_dmr_test
    LIBS parameter base device hcontainer
    SOURCES ../module_dm/test/test_write_dmr.cpp
            ../module_dm/write_dmr.cpp
            ../../source_cell/ucell_io.cpp
            ../../source_basis/module_ao/parallel_orbitals.cpp
            ../../source_hamilt/module_hcontainer/test/tmp_mocks.cpp)
endif()

addtest(
  TARGET MODULE_IO_read_wfc_nao_test
  LIBS parameter base device
  SOURCES test_read_wfc_nao.cpp ../module_wf/read_wfc_nao.cpp
          ../module_wf/write_wfc_nao.cpp ../../source_psi/psi.cpp
          ../../source_basis/module_ao/parallel_orbitals.cpp)

add_test(
  NAME MODULE_IO_read_wfc_nao_test_parallel
  COMMAND mpirun -np 2 "$<TARGET_FILE:MODULE_IO_read_wfc_nao_test>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_LIBRI)
  addtest(
    TARGET MODULE_IO_restart_exx_csr_test
    LIBS parameter base device
    SOURCES test_restart_exx_csr.cpp ../module_restart/restart_exx_csr.cpp
            ../module_hs/hs_sparse_io.cpp ../module_hs/lat_r_csr.cpp
            tmp_mocks.cpp)
endif()

addtest(TARGET MODULE_IO_write_elf_logic_test SOURCES test_write_elf_logic.cpp)
