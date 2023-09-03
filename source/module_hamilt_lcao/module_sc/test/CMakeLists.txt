if(ENABLE_LCAO)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
AddTest(
  TARGET sc_test
  LIBS ${math_libs} base device
  SOURCES sc_test.cpp ../spin_constrain.cpp ../sc_parse_json.cpp
)
endif() 
