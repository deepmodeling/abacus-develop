remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)
remove_definitions(-D__LCAO)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES unitcell_test_pw_para.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_CELL_unitcell_test_pw
  LIBS parameter ${math_libs} base device
  SOURCES unitcell_test_pw.cpp ../unitcell.cpp ../read_atoms.cpp ../atom_spec.cpp ../update_cell.cpp ../bcast_cell.cpp
	../atom_pseudo.cpp ../pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp
  ../read_stru.cpp   ../read_atom_species.cpp
	../read_pp_vwr.cpp ../read_pp_blps.cpp ../../source_io/output.cpp
  ../../source_estate/read_pseudo.cpp ../../source_estate/cal_nelec_nband.cpp
  ../../source_estate/read_orb.cpp    ../print_cell.cpp
  ../../source_estate/cal_wfc.cpp  ../sep.cpp   ../sep_cell.cpp
)

find_program(BASH bash)
add_test(NAME MODULE_CELL_unitcell_test_pw_para
      COMMAND ${BASH} unitcell_test_pw_para.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
