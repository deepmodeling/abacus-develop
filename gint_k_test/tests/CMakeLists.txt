# Add Google Test
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Add unit tests
add_executable(gint_k_test gint_k_test.cpp)
target_link_libraries(gint_k_test ${GTEST_LIBRARIES} pthread gint_k_lib)

# Add benchmark tests
add_executable(benchmark_gint_k ../benchmark_gint_k.cpp)
target_link_libraries(benchmark_gint_k benchmark::benchmark gint_k_lib)

