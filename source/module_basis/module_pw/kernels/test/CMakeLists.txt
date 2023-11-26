add_definitions(-D__NORMAL)

AddTest(
    TARGET PW_Kernels_UTs
    LIBS ${math_libs} psi device
    SOURCES pw_op_test.cpp
    ../../../../module_base/tool_quit.cpp ../../../../module_base/global_variable.cpp
    ../../../../module_base/para_json.cpp
    ../../../../module_base/parallel_global.cpp ../../../../module_base/parallel_reduce.cpp
    ../../../../module_base/complexmatrix.cpp ../../../../module_base/matrix.cpp ../../../../module_base/memory.cpp
    ../../../../module_base/libm/branred.cpp ../../../../module_base/libm/sincos.cpp
)