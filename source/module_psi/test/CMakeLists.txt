AddTest(
    TARGET psi_UT
    LIBS ${math_libs} base device
    SOURCES 
        psi_test.cpp
        ../psi.cpp 
)

if(ENABLE_LCAO)
AddTest(
    TARGET psi_initializer_unit_test
    LIBS ${math_libs} base device psi psi_initializer planewave
    SOURCES
        psi_initializer_unit_test.cpp
        ../../module_hamilt_pw/hamilt_pwdft/soc.cpp
        ../../module_cell/atom_spec.cpp
        ../../module_cell/parallel_kpoints.cpp
        mock_unitcell.cpp
)
endif()

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})