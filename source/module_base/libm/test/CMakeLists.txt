remove_definitions(-D__MPI)

AddTest(
    TARGET Base_Libm_UTs
    SOURCES libm_test.cpp
    ../branred.cpp
    ../cexp.cpp
    ../exp.cpp
    ../sincos.cpp
)

