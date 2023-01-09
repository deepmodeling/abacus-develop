AddTest(
    TARGET Hamilt_Kernels_UTs
    LIBS ${math_libs} device base
    SOURCES ekinetic_op_test.cpp nonlocal_op_test.cpp veff_op_test.cpp meta_op_test.cpp
    ../../../src_parallel/parallel_reduce.cpp
    ../../../src_parallel/parallel_global.cpp
)
