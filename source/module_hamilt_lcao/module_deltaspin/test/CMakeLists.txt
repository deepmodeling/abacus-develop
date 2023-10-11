if(ENABLE_LCAO)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET deltaspin_parse_json_test
  LIBS ${math_libs} base device
  SOURCES parse_json_test.cpp ../spin_constrain.cpp ../sc_parse_json.cpp
)

AddTest(
  TARGET deltaspin_spin_constrain_test
  LIBS ${math_libs} base device
  SOURCES spin_constrain_test.cpp ../spin_constrain.cpp ../sc_parse_json.cpp
)

AddTest(
  TARGET deltaspin_basic_func_test
  LIBS ${math_libs} base device
  SOURCES basic_test.cpp ../basic_funcs.cpp
)
endif() 
