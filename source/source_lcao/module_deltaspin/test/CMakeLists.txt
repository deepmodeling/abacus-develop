remove_definitions(-D__CUDA)

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
  LIBS ${math_libs} base
  SOURCES deltaspin_core_test.cpp
)

AddTest(
  TARGET deltaspin_lambda_update_strategies_test
  LIBS ${math_libs} base device parameter
  SOURCES lambda_update_strategies_test.cpp
    ../lambda_update_strategies.cpp
)

# TODO: lambda_loop_helper_test.cpp needs API update to match current SpinConstrain<TK> signature
# AddTest(
#   TARGET deltaspin_lambda_loop_helper_test
#   LIBS ${math_libs} base device parameter
#   SOURCES lambda_loop_helper_test.cpp
#     ../spin_constrain.cpp
#     ../lambda_loop_helper.cpp
# )

endif()
