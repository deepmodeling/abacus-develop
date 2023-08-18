
if(ENABLE_LCAO)
AddTest(
  TARGET deltaspin_basic
  LIBS ${math_libs} base device
  SOURCES basic_test.cpp ../basic_funcs.cpp
)
endif() 
