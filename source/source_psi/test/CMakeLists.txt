AddTest(
    TARGET MODULE_PSI_Unittests
    LIBS parameter ${math_libs} base device 
    SOURCES 
        psi_test.cpp
        ../psi.cpp 
)

if(ENABLE_LCAO)
AddTest(
    TARGET MODULE_PSI_initializer_unit_test
    LIBS parameter ${math_libs} base device psi psi_initializer planewave 
    SOURCES
        psi_initializer_unit_test.cpp
        ../../source_pw/module_pwdft/soc.cpp
        ../../source_cell/atom_spec.cpp
        ../../source_cell/parallel_kpoints.cpp
        ../../source_cell/test/support/mock_unitcell.cpp
        ../../source_io/module_output/orb_io.cpp
        ../../source_io/module_output/write_pao.cpp
        ../../source_io/module_wf/read_wfc_pw.cpp
        ../../source_io/module_output/binstream.cpp
        ../../source_io/module_output/filename.cpp
)
endif()

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
