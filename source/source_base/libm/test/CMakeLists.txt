abacus_disable_feature_definitions(__MPI)

AddTest(
    TARGET MODULE_BASE_LIBM_UTs
    SOURCES libm_test.cpp
    ../branred.cpp
    ../cexp.cpp
    ../exp.cpp
    ../sincos.cpp
)

