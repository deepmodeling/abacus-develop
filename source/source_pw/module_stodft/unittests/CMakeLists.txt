abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

# Regression test for the occupation / entropy kernels of the SDFT solver:
# `1 - f` must not be computed by subtraction and the expanded function must
# not contain an artificial cutoff (|x| > 36 used to return exactly 0).
AddTest(TARGET MODULE_PW_STO_FUNC_UTs LIBS parameter base device SOURCES
        test_sto_func.cpp ../sto_func.cpp ../../../source_estate/occupy.cpp)

# Regression test for the pure DOS quality guards: the Chebyshev-window
# sanitiser and the DOS parameter / result checks.
AddTest(TARGET MODULE_PW_STO_DOS_QUALITY_UTs LIBS base device SOURCES
        test_sto_dos_quality.cpp ../sto_dos_quality.cpp)
