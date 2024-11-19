if(ENABLE_LCAO)

AddTest(
  TARGET deltaspin_basic_func_test
  LIBS ${math_libs} base device parameter
  SOURCES basic_test.cpp
    ../basic_funcs.cpp
)

AddTest(
  TARGET deltaspin_spin_constrain_test
  LIBS ${math_libs} base device parameter
  SOURCES spin_constrain_test.cpp
    ../spin_constrain.cpp
    ../../../module_basis/module_ao/parallel_2d.cpp
    ../../../module_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET deltaspin_template_helpers
  LIBS ${math_libs} base device parameter
  SOURCES template_helpers_test.cpp
    ../spin_constrain.cpp
    ../template_helpers.cpp
)
endif() 
