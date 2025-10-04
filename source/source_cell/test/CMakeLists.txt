remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

find_program(BASH bash)
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
    ../../source_estate/read_pseudo.cpp
    ../../source_estate/cal_wfc.cpp
    ../../source_estate/cal_nelec_nband.cpp
    ../../source_estate/read_orb.cpp
    ../sep.cpp
    ../sep_cell.cpp
)

add_library(cell_info OBJECT ${cell_simple_srcs})

AddTest(
  TARGET MODULE_CELL_read_pp
  LIBS parameter  ${math_libs} base device
  SOURCES read_pp_test.cpp ../atom_pseudo.cpp ../pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp ../read_pp_vwr.cpp ../read_pp_blps.cpp
  ../../source_io/output.cpp
)

AddTest(
  TARGET MODULE_CELL_pseudo_nc
  LIBS parameter  ${math_libs} base device
  SOURCES pseudo_nc_test.cpp ../pseudo.cpp ../atom_pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp ../read_pp_vwr.cpp
  	../read_pp_blps.cpp ../../source_io/output.cpp
)

AddTest(
  TARGET MODULE_CELL_atom_pseudo
  LIBS parameter  ${math_libs} base device
  SOURCES atom_pseudo_test.cpp ../atom_pseudo.cpp ../pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp
  	../read_pp_vwr.cpp ../read_pp_blps.cpp ../../source_io/output.cpp
)

AddTest(
  TARGET MODULE_CELL_atom_spec
  LIBS parameter ${math_libs} base device
  SOURCES atom_spec_test.cpp ../atom_spec.cpp ../atom_pseudo.cpp ../pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp
  	../read_pp_upf100.cpp ../read_pp_vwr.cpp ../read_pp_blps.cpp ../../source_io/output.cpp
)

AddTest(
  TARGET MODULE_CELL_klist_test
  LIBS parameter  ${math_libs} base device symmetry
  SOURCES klist_test.cpp ../klist.cpp ../parallel_kpoints.cpp ../../source_io/output.cpp ../k_vector_utils.cpp
)

AddTest(
  TARGET MODULE_CELL_klist_test_para1
  LIBS parameter  ${math_libs} base device symmetry
  SOURCES klist_test_para.cpp ../klist.cpp ../parallel_kpoints.cpp ../../source_io/output.cpp ../k_vector_utils.cpp
)

add_test(NAME MODULE_CELL_klist_test_para4
      COMMAND ${BASH} klist_test_para.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_CELL_ParaKpoints
  LIBS parameter  MPI::MPI_CXX
  SOURCES parallel_kpoints_test.cpp ../../source_base/global_variable.cpp ../../source_base/parallel_global.cpp
  	../../source_base/parallel_common.cpp ../../source_base/parallel_comm.cpp ../parallel_kpoints.cpp
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
  LIBS parameter  ${math_libs} base device cell_info  symmetry
  SOURCES unitcell_test.cpp ../../source_io/output.cpp ../../source_estate/cal_ux.cpp

)

AddTest(
  TARGET MODULE_CELL_unitcell_test_readpp
  LIBS parameter  ${math_libs} base device cell_info
  SOURCES unitcell_test_readpp.cpp ../../source_io/output.cpp
)

AddTest(
  TARGET MODULE_CELL_unitcell_test_para
  LIBS parameter  ${math_libs} base device cell_info
  SOURCES unitcell_test_para.cpp ../../source_io/output.cpp
)

AddTest(
  TARGET MODULE_CELL_unitcell_test_setupcell
  LIBS parameter  ${math_libs} base device cell_info
  SOURCES unitcell_test_setupcell.cpp ../../source_io/output.cpp
)

add_test(NAME MODULE_CELL_unitcell_test_parallel
      COMMAND ${BASH} unitcell_test_parallel.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_CELL_index_test
  LIBS parameter  ${math_libs} base device
  SOURCES cell_index_test.cpp ../cell_index.cpp
)

AddTest(
  TARGET MODULE_CELL_SEP_TEST
  LIBS parameter  ${math_libs} base device
  SOURCES read_sep_test.cpp ../sep.cpp
)

add_test(NAME MODULE_CELL_read_sep_parallel
      COMMAND ${BASH} bcast_read_sep_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_CELL_SEP_CELL_TEST
  LIBS parameter  ${math_libs} base device
  SOURCES sepcell_test.cpp ../sep.cpp ../sep_cell.cpp
)

add_test(NAME MODULE_CELL_sep_cell_parallel
      COMMAND ${BASH} bcast_sep_cell_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
