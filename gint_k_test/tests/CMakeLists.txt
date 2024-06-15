# 添加 Google Test
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# 添加单元测试
add_executable(gint_k_test gint_k_test.cpp)
target_link_libraries(gint_k_test ${GTEST_LIBRARIES} pthread gint_k_lib)

# 添加性能测试
add_executable(benchmark_gint_k ../benchmark_gint_k.cpp)
target_link_libraries(benchmark_gint_k benchmark::benchmark gint_k_lib)

