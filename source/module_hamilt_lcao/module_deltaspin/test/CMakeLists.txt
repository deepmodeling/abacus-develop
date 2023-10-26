if(ENABLE_LCAO)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET deltaspin_parse_json_test
  LIBS ${math_libs} base device
  SOURCES parse_json_test.cpp
    ../spin_constrain.cpp
    ../sc_parse_json.cpp
)

AddTest(
  TARGET deltaspin_basic_func_test
  LIBS ${math_libs} base device
  SOURCES basic_test.cpp
    ../basic_funcs.cpp
)

AddTest(
  TARGET deltaspin_spin_constrain_test
  LIBS ${math_libs} base device
  SOURCES spin_constrain_test.cpp
    ../spin_constrain.cpp
    ../sc_parse_json.cpp
    ../../../module_basis/module_ao/parallel_2d.cpp
    ../../../module_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET deltaspin_template_helpers
  LIBS ${math_libs} base device
  SOURCES template_helpers_test.cpp
    ../spin_constrain.cpp
    ../template_helpers.cpp
)

AddTest(
  TARGET deltaspin_cal_h_lambda
  LIBS ${math_libs} base device
  SOURCES cal_h_lambda_test.cpp
    ../spin_constrain.cpp
    ../cal_h_lambda.cpp
    ../../../module_basis/module_ao/parallel_2d.cpp
    ../../../module_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET deltaspin_bcast_scdata
  LIBS ${math_libs} base device
  SOURCES bcast_scdata_test.cpp
    ../spin_constrain.cpp
    ../sc_parse_json.cpp
)

install(FILES bcast_scdata_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(BASH bash)
add_test(NAME deltaspin_bcast_scdata_para
      COMMAND ${BASH} bcast_scdata_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET deltaspin_init_sc
  LIBS ${math_libs} base device cell_info
  SOURCES init_sc_test.cpp
    ../spin_constrain.cpp
    ../sc_parse_json.cpp
    ../init_sc.cpp
    ../../../module_io/output.cpp
    ../../../module_basis/module_ao/parallel_2d.cpp
    ../../../module_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET deltaspin_cal_mw_helper
  LIBS ${math_libs} base device
  SOURCES cal_mw_helper_test.cpp
    ../spin_constrain.cpp
    ../cal_mw_helper.cpp
)
endif() 
