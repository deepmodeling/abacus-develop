cmake_minimum_required(VERSION 3.10)
project(UpsiTest)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(UpsiTest upsi_test.cpp)
target_include_directories(UpsiTest PRIVATE ../../../module_basis/module_ao)
target_link_libraries(UpsiTest ${GTEST_BOTH_LIBRARIES} pthread)

