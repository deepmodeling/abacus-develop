AddTest(
    TARGET XC_Functional_UTs
    LIBS parameter ${math_libs} device base container 
    SOURCES xc_functional_op_test.cpp
)
