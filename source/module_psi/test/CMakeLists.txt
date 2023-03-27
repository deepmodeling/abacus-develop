AddTest(
    TARGET psi_UT
    LIBS ${math_libs} base device
    SOURCES 
        psi_test.cpp
        ../psi.cpp 
)