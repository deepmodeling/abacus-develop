abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__LCAO)

file(COPY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
configure_file(unitcell_test_pw_para.sh ${CMAKE_CURRENT_BINARY_DIR}/unitcell_test_pw_para.sh COPYONLY)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES unitcell_test_pw_para.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_CELL_unitcell_test_pw
  LIBS parameter base device
  SOURCES unitcell_test_pw.cpp ../unitcell.cpp ../read_atoms.cpp ../read_atoms_helper.cpp ../atom_spec.cpp ../update_cell.cpp ../bcast_cell.cpp
	../atom_pseudo.cpp ../pseudo.cpp ../read_pp.cpp ../read_pp_complete.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp
  ../read_stru.cpp   ../read_atom_species.cpp
	../read_pp_vwr.cpp ../read_pp_blps.cpp
  ../../source_estate/read_pseudo.cpp ../../source_estate/cal_nelec_nband.cpp
  ../../source_estate/read_orb.cpp    ../print_cell.cpp
  ../../source_estate/cal_wfc.cpp  ../sep.cpp   ../sep_cell.cpp
)

find_program(BASH bash)
add_test(NAME MODULE_CELL_unitcell_test_pw_para
      COMMAND ${BASH} unitcell_test_pw_para.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
