abacus_disable_feature_definitions(__CUDA)

if(ENABLE_LCAO)

AddTest(
  TARGET MODULE_LCAO_deltaspin_basic_func_test
  LIBS ${math_libs} base device parameter
  SOURCES basic_test.cpp
    ../basic_funcs.cpp
)

AddTest(
  TARGET MODULE_LCAO_deltaspin_spin_constrain_test
  LIBS ${math_libs} base device parameter
  SOURCES spin_constrain_test.cpp
    ../spin_constrain.cpp
    ../../../source_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET MODULE_LCAO_deltaspin_template_helpers
  LIBS ${math_libs} base device parameter
  SOURCES template_helpers_test.cpp
    ../spin_constrain.cpp
    ../template_helpers.cpp
)

AddTest(
  TARGET deltaspin_pw_test
  LIBS ${math_libs} base device parameter
  SOURCES deltaspin_pw_test.cpp
)

AddTest(
  TARGET deltaspin_core_test
  LIBS ${math_libs} base device
  SOURCES deltaspin_core_test.cpp
)

endif()
