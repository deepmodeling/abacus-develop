AddTest(
    TARGET deltaspin_force_cpu_test
    SOURCES test_force_op.cpp ../force_op.cpp
)

if(USE_CUDA)
AddTest(
    TARGET deltaspin_force_gpu_test
    LIBS CUDA::cudart
    KEEP_FEATURE_DEFINITIONS __CUDA
    SOURCES test_force_op.cpp ../force_op.cpp ../cuda/force_op.cu
)
endif()

abacus_disable_feature_definitions(__LCAO)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)

AddTest(
    TARGET MODULE_PW_Hamilt_Kernels_UTs
    LIBS parameter device base container
    SOURCES ekinetic_op_test.cpp nonlocal_op_test.cpp veff_op_test.cpp meta_op_test.cpp
    force_op_test.cpp stress_op_test.cpp wf_op_test.cpp vnl_op_test.cpp
    stress_op_mgga_test.cpp
)

AddTest(
    TARGET onsite_op_test
    LIBS parameter device base container
    SOURCES onsite_op_test.cpp
)
