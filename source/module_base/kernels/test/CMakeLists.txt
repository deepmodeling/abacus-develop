remove_definitions(-D__MPI)

AddTest(
    TARGET MODULE_BASE_KERNELS_Unittests
    LIBS parameter ${math_libs} base device 
    SOURCES math_ylm_op_test.cpp math_kernel_test.cpp
)
