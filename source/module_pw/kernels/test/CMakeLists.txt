add_definitions(-D__NORMAL)

AddTest(
    TARGET PW_Kernels_UTs
    LIBS ${math_libs} psi device
    SOURCES pw_op_test.cpp
    ../../../module_base/tool_quit.cpp ../../../module_base/global_variable.cpp
    ../../../src_parallel/parallel_global.cpp ../../../src_parallel/parallel_reduce.cpp
    ../../../module_base/complexmatrix.cpp ../../../module_base/matrix.cpp
)