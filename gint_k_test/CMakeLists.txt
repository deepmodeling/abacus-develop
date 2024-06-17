cmake_minimum_required(VERSION 3.10)
project(GintKTest)

# Add Google Benchmark and Google Test
add_subdirectory(../benchmark benchmark_build)
include_directories(${benchmark_SOURCE_DIR}/include)
add_subdirectory(tests)

# Add files
add_library(gint_k_lib gint_k.cpp)

