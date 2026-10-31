enable_testing()

# TODO: Enable all tests after DeltaSpin + DFT+U code is merged into develop.
# These tests require the DeltaSpin module (source/module_hamilt_general/hamilt_deltaspin.h)
# and associated DFT+U with spin constraint functionality.
# Uncomment the add_test blocks below once the feature branch is merged.

find_program(BASH bash)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(
  NAME 17_DS_DFTU_onsite_mesh
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_onsite_proj_mesh.py
          ${ABACUS_BIN_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/50_FeO_O_first_Fe_second)
if(ENABLE_ASAN)
  add_test(
    NAME 17_DS_DFTU_test_with_asan
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/17_DS_DFTU)
else()
  add_test(
    NAME 17_DS_DFTU
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/17_DS_DFTU)
endif()
