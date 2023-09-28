AddTest(
    TARGET psi_UT
    LIBS ${math_libs} base device
    SOURCES 
        psi_test.cpp
        ../psi.cpp 
)

if(ENABLE_LCAO)
AddTest(
    TARGET psi_initializer_integrate_test
    SOURCES
        psi_initializer_integrated_test.cpp
)
endif()

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})