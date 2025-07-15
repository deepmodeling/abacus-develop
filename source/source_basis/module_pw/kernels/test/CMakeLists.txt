add_definitions(-D__NORMAL)

AddTest(
    TARGET PW_Kernels_UTs
    LIBS parameter ${math_libs} psi device
    SOURCES pw_op_test.cpp
    ../../../../source_base/tool_quit.cpp ../../../../source_base/global_variable.cpp
    ../../../../source_base/parallel_global.cpp ../../../../source_base/parallel_reduce.cpp
    ../../../../source_base/parallel_comm.cpp
    ../../../../source_base/complexmatrix.cpp ../../../../source_base/matrix.cpp ../../../../source_base/memory.cpp
    ../../../../source_base/libm/branred.cpp ../../../../source_base/libm/sincos.cpp
    ../../../../source_base/module_external/blas_connector_base.cpp ../../../../source_base/module_external/blas_connector_vector.cpp ../../../../source_base/module_external/blas_connector_matrix.cpp
)