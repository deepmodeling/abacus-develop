remove_definitions(-D__MPI)

AddTest(
    TARGET Base_Kernels_UTs
    LIBS ${math_libs} base device
    SOURCES math_op_test.cpp
)

