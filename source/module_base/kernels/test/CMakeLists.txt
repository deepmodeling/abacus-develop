remove_definitions(-D__MPI)

AddTest(
    TARGET Base_Kernels_UTs
    LIBS ${math_libs} base device
    SOURCES math_op_test.cpp ../../../src_parallel/parallel_reduce.cpp
    ../../../src_parallel/parallel_kpoints.cpp
    ../../../src_parallel/parallel_global.cpp
    ../../../src_parallel/parallel_common.cpp
)

