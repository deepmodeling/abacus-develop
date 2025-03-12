remove_definitions(-D__MPI)

AddTest(
    TARGET Base_Kernels_UTs
    LIBS parameter ${math_libs} base device 
    SOURCES math_op_test.cpp math_kernel_test.cpp
)
