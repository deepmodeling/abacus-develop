find_program(BASH bash)
find_package(Python3 REQUIRED)

add_executable(sum_cube ${ABACUS_TEST_DIR}/integrate/tools/sum_cube.cpp)
install(TARGETS sum_cube DESTINATION ${ABACUS_TEST_DIR}/integrate/tools/)
message(STATUS "Target sum_cube has been defined in ${CMAKE_CURRENT_SOURCE_DIR}") # 添加此行
    
if(ENABLE_COVERAGE)
  add_coverage(sum_cube)
endif()

if(ENABLE_ASAN)
    add_test(
        NAME integrated_test_with_asan
        COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
    )

else()
    add_test(
        NAME integrated_test
        COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
    )
endif()