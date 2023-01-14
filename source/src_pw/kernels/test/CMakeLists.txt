remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

AddTest(
    TARGET SRC_PW_Kernels_UTs
    LIBS ${math_libs} base device
    SOURCES force_op_test.cpp stress_op_test.cpp wf_op_test.cpp vnl_op_test.cpp
    ../../../src_parallel/parallel_reduce.cpp
    ../../../src_parallel/parallel_global.cpp
)