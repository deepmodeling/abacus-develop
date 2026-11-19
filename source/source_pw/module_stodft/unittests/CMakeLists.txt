# Sto_Func is a scalar CPU utility; its tests also run in serial builds.
abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

addtest(TARGET MODULE_PW_sto_func SOURCES test_sto_func.cpp ../sto_func.cpp)
