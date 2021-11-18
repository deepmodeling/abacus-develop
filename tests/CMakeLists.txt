add_subdirectory(integrate)

set(CMAKE_LINK_WHAT_YOU_USE TRUE)
find_package(Threads)

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG "origin/main"
)
FetchContent_MakeAvailable(googletest)
include_directories(${ABACUS_SOURCE_DIR})
function(AddTest)
  cmake_parse_arguments(UT "DYN" "TARGET" "LIBS;DYN_LIBS;STATIC_LIBS;SOURCES;DEPENDS" ${ARGN})
  add_executable(${UT_TARGET} ${UT_SOURCES})
  #dependencies & link library
  target_link_libraries(${UT_TARGET} ${UT_LIBS} ${ABACUS_LINK_LIBRARIES} pthread gtest_main)
  install(TARGETS ${UT_TARGET} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests )
  add_test(NAME ${UT_TARGET}
  COMMAND ${UT_TARGET}
  WORKING_DIRECTORY $<TARGET_FILE_DIR:${UT_TARGET}>
  )
endfunction(AddTest)

get_target_property(ABACUS_LINK_LIBRARIES ${ABACUS_BIN_NAME} LINK_LIBRARIES)

set(TEST_BASE_DIR ${CMAKE_SOURCE_DIR}/source/module_base)
include_directories(${TEST_BASE_DIR})

AddTest(
  TARGET base_matrix3
  SOURCES ${ABACUS_TEST_DIR}/module_base/test_matrix3.cpp
)
AddTest(
  TARGET base_blas_connector
  SOURCES ${ABACUS_TEST_DIR}/module_base/test_blas_connector.cpp
)
