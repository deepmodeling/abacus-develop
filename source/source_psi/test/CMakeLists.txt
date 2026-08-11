AddTest(
    TARGET MODULE_PSI_Unittests
    LIBS parameter base device
    SOURCES 
        psi_test.cpp
        ../psi.cpp 
)

if(ENABLE_LCAO)
AddTest(
    TARGET MODULE_PSI_init_test
    LIBS parameter base device psi psi_init planewave
    SOURCES
        psi_init_test.cpp
        ../../source_pw/module_pwdft/soc.cpp
        ../../source_cell/atom_spec.cpp
        ../../source_cell/test/support/mock_unitcell.cpp
        
        ../../source_io/module_output/write_pao.cpp
        ../../source_io/module_wf/read_wfc_pw.cpp
)
endif()

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
