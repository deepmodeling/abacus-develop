remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)

AddTest(
    TARGET Hamilt_Kernels_UTs
    LIBS ${math_libs} device base container
    SOURCES ekinetic_op_test.cpp nonlocal_op_test.cpp veff_op_test.cpp meta_op_test.cpp
    force_op_test.cpp stress_op_test.cpp wf_op_test.cpp vnl_op_test.cpp
    stress_op_mgga_test.cpp
)
