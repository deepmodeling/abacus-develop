remove_definitions(-D__MPI)

AddTest(
    TARGET MODULE_BASE_LIBM_UTs
    SOURCES libm_test.cpp
    ../branred.cpp
    ../cexp.cpp
    ../exp.cpp
    ../sincos.cpp
)

