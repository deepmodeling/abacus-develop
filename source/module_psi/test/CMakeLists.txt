AddTest(
    TARGET psi_UT
    LIBS parameter ${math_libs} base device 
    SOURCES 
        psi_test.cpp
        ../psi.cpp 
)

if(ENABLE_LCAO)
AddTest(
    TARGET psi_initializer_unit_test
    LIBS parameter ${math_libs} base device psi psi_initializer planewave 
    SOURCES
        psi_initializer_unit_test.cpp
        ../../module_hamilt_pw/hamilt_pwdft/soc.cpp
        ../../module_cell/atom_spec.cpp
        ../../module_cell/parallel_kpoints.cpp
        ../../module_cell/test/support/mock_unitcell.cpp
        ../../module_io/orb_io.cpp
        ../../module_io/write_pao.cpp
        ../../module_io/read_wfc_pw.cpp
        ../../module_io/binstream.cpp
)
endif()

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})