abacus_disable_feature_definitions(__LCAO)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)

AddTest(
    TARGET MODULE_PW_Hamilt_Kernels_UTs
    LIBS parameter ${math_libs} device base container
    SOURCES ekinetic_op_test.cpp nonlocal_op_test.cpp veff_op_test.cpp meta_op_test.cpp
    force_op_test.cpp stress_op_test.cpp wf_op_test.cpp vnl_op_test.cpp
    stress_op_mgga_test.cpp
)

AddTest(
    TARGET onsite_op_test
    LIBS parameter ${math_libs} device base container
    SOURCES onsite_op_test.cpp
)
