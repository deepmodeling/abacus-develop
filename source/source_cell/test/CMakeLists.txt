abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

find_program(BASH bash)
file(COPY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
configure_file(bcast_atom_pseudo_test.sh ${CMAKE_CURRENT_BINARY_DIR}/bcast_atom_pseudo_test.sh COPYONLY)
configure_file(bcast_atom_spec_test.sh ${CMAKE_CURRENT_BINARY_DIR}/bcast_atom_spec_test.sh COPYONLY)
configure_file(parallel_kpoints_test.sh ${CMAKE_CURRENT_BINARY_DIR}/parallel_kpoints_test.sh COPYONLY)
configure_file(klist_test_para.sh ${CMAKE_CURRENT_BINARY_DIR}/klist_test_para.sh COPYONLY)
configure_file(unitcell_test_parallel.sh ${CMAKE_CURRENT_BINARY_DIR}/unitcell_test_parallel.sh COPYONLY)
configure_file(bcast_read_sep_test.sh ${CMAKE_CURRENT_BINARY_DIR}/bcast_read_sep_test.sh COPYONLY)
configure_file(bcast_sep_cell_test.sh ${CMAKE_CURRENT_BINARY_DIR}/bcast_sep_cell_test.sh COPYONLY)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES bcast_atom_pseudo_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES bcast_atom_spec_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES parallel_kpoints_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES klist_test_para.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES unitcell_test_parallel.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES bcast_read_sep_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES bcast_sep_cell_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND cell_simple_srcs
    ../unitcell.cpp
    ../update_cell.cpp
    ../bcast_cell.cpp
    ../print_cell.cpp
    ../read_stru.cpp
    ../read_atom_species.cpp
    ../read_atoms.cpp
    ../read_atoms_helper.cpp
    ../atom_spec.cpp
    ../atom_pseudo.cpp
    ../pseudo.cpp
    ../read_pp.cpp
    ../read_pp_complete.cpp
    ../read_pp_upf201.cpp
    ../read_pp_upf100.cpp
    ../read_pp_vwr.cpp
    ../read_pp_blps.cpp
    ../check_atomic_stru.cpp
    ../read_pp_ucell.cpp
    ../cal_wfc.cpp
    ../cal_nelec_nband.cpp
    ../read_orb.cpp
    ../sep.cpp
    ../sep_cell.cpp
    ../cell_tools.cpp
)

add_library(cell_info OBJECT ${cell_simple_srcs})

AddTest(
  TARGET MODULE_CELL_read_pp
  LIBS base device
  SOURCES read_pp_test.cpp ../atom_pseudo.cpp ../pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp ../read_pp_vwr.cpp ../read_pp_blps.cpp

)

AddTest(
  TARGET MODULE_CELL_pseudo_nc
  LIBS base device
  SOURCES pseudo_nc_test.cpp ../pseudo.cpp ../atom_pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp ../read_pp_vwr.cpp
  	../read_pp_blps.cpp
)

AddTest(
  TARGET MODULE_CELL_atom_pseudo
  LIBS base device
  SOURCES atom_pseudo_test.cpp ../atom_pseudo.cpp ../pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp
  	../read_pp_vwr.cpp ../read_pp_blps.cpp
)

AddTest(
  TARGET MODULE_CELL_atom_spec
  LIBS base device
  SOURCES atom_spec_test.cpp ../atom_spec.cpp ../atom_pseudo.cpp ../pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp
  	../read_pp_upf100.cpp ../read_pp_vwr.cpp ../read_pp_blps.cpp
)

AddTest(
  TARGET MODULE_CELL_klist_test
  LIBS base device symmetry
  SOURCES klist_test.cpp ../klist.cpp ../parallel_kpoints.cpp ../k_vector_utils.cpp ../reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_CELL_reciprocal_grid_test
  LIBS base device symmetry
  SOURCES reciprocal_grid_test.cpp ../reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_CELL_qlist_test
  LIBS base device symmetry
  SOURCES qlist_test.cpp ../qlist.cpp ../reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_CELL_little_group_test
  LIBS base device symmetry
  SOURCES little_group_test.cpp
)

AddTest(
  TARGET MODULE_CELL_klist_test_para1
  LIBS base device symmetry
  SOURCES klist_test_para.cpp ../klist.cpp ../parallel_kpoints.cpp ../k_vector_utils.cpp ../reciprocal_grid.cpp
)

add_test(NAME MODULE_CELL_klist_test_para4
      COMMAND ${BASH} klist_test_para.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_CELL_ParaKpoints
  LIBS MPI::MPI_CXX
  SOURCES parallel_kpoints_test.cpp ../../source_base/global_variable.cpp ../../source_base/parallel_global.cpp
  	../../source_base/parallel_common.cpp ../../source_base/parallel_comm.cpp ../parallel_kpoints.cpp ../../source_base/tool_quit.cpp ../../source_base/global_variable.cpp ../../source_base/global_file.cpp ../../source_base/global_function.cpp ../../source_base/memory_recorder.cpp ../../source_base/timer.cpp ../../source_base/parallel_reduce.cpp
)

# Add unit test for read_atoms_helper
AddTest(
  TARGET MODULE_CELL_read_atoms_helper_test
  LIBS base device
  SOURCES read_atoms_hlp_test.cpp
    ../read_atoms.cpp
    ../read_atoms_helper.cpp
    ../cell_tools.cpp
    ../read_orb.cpp
    ../read_stru.cpp
    ../print_cell.cpp
    ../atom_spec.cpp
    ../unitcell.cpp
    ../update_cell.cpp
    ../bcast_cell.cpp
    ../atom_pseudo.cpp
    ../pseudo.cpp
    ../read_pp.cpp
    ../read_pp_complete.cpp
    ../read_pp_upf201.cpp
    ../read_pp_upf100.cpp
    ../read_pp_vwr.cpp
    ../read_pp_blps.cpp
    ../read_atom_species.cpp
    ../sep.cpp
    ../sep_cell.cpp
)

find_program(BASH bash)

add_test(NAME MODULE_CELL_bcast_atom_pseudo_test
      COMMAND ${BASH} bcast_atom_pseudo_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME MODULE_CELL_bcast_atom_spec_test
      COMMAND ${BASH} bcast_atom_spec_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME MODULE_CELL_parallel_kpoints_test
      COMMAND ${BASH} parallel_kpoints_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_CELL_unitcell_test
  LIBS base device cell_info symmetry
  SOURCES unitcell_test.cpp ../cal_ux.cpp

)

AddTest(
  TARGET MODULE_CELL_unitcell_test_readpp
  LIBS base device cell_info
  SOURCES unitcell_test_readpp.cpp
)

AddTest(
  TARGET MODULE_CELL_unitcell_test_para
  LIBS base device cell_info
  SOURCES unitcell_test_para.cpp
)

AddTest(
  TARGET MODULE_CELL_unitcell_test_setupcell
  LIBS base device cell_info
  SOURCES unitcell_test_setupcell.cpp
)

add_test(NAME MODULE_CELL_unitcell_test_parallel
      COMMAND ${BASH} unitcell_test_parallel.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_CELL_index_test
  LIBS base device
  SOURCES cell_index_test.cpp ../cell_index.cpp
)

AddTest(
  TARGET MODULE_CELL_SEP_TEST
  LIBS base device
  SOURCES read_sep_test.cpp ../sep.cpp
)

add_test(NAME MODULE_CELL_read_sep_parallel
      COMMAND ${BASH} bcast_read_sep_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_CELL_SEP_CELL_TEST
  LIBS base device
  SOURCES sepcell_test.cpp ../sep.cpp ../sep_cell.cpp
)

add_test(NAME MODULE_CELL_sep_cell_parallel
      COMMAND ${BASH} bcast_sep_cell_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_CELL_magnetism
  LIBS base device
  SOURCES magnetism_test.cpp ../magnetism.cpp
)

AddTest(
  TARGET MODULE_CELL_cif_io
  LIBS parameter base device
  SOURCES cif_io_test.cpp ../cif_io.cpp
)

add_test(NAME MODULE_CELL_cif_io_parallel
      COMMAND mpirun -np 4 ./MODULE_CELL_cif_io
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
