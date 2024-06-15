cmake_minimum_required(VERSION 3.10)
project(GintKTest)

# 添加 Google Benchmark 和 Google Test
add_subdirectory(../benchmark benchmark_build)
include_directories(${benchmark_SOURCE_DIR}/include)
add_subdirectory(tests)

# 添加主程序源文件
add_library(gint_k_lib gint_k.cpp)

