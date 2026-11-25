addtest(
  TARGET MODULE_PSI_Unittests
  LIBS parameter base device
  SOURCES psi_test.cpp ../psi.cpp)

if(ENABLE_LCAO)
  addtest(
    TARGET MODULE_PSI_init_test
    LIBS parameter base device psi psi_init planewave
    SOURCES psi_init_unit_test.cpp
            ../../source_cell/soc.cpp
            ../../source_cell/atom_spec.cpp
            ../../source_cell/test/support/mock_unitcell.cpp
            ../../source_io/module_output/write_pao.cpp
            ../../source_io/module_wf/read_wfc_pw.cpp)
endif()

abacus_copy_test_files(support)
