cmake_minimum_required(VERSION 3.16)
project(nep_benchmark CXX)

option(USE_OPENMP "Enable OpenMP" ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_path(NEP_INCLUDE_DIR nep.h
  HINTS ${NEP_DIR}
  PATH_SUFFIXES include src
  NO_DEFAULT_PATH)
find_library(NEP_LIBRARY nep
  HINTS ${NEP_DIR}
  PATH_SUFFIXES lib .
  NO_DEFAULT_PATH)

if(NOT NEP_INCLUDE_DIR)
  set(NEP_INCLUDE_DIR "${NEP_DIR}/../src")
  if(NOT EXISTS "${NEP_INCLUDE_DIR}/nep.h")
    set(NEP_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../third_party/NEP_CPU/src")
  endif()
endif()

if(NOT NEP_LIBRARY)
  set(NEP_LIBRARY "${NEP_DIR}/libnep.so")
  if(NOT EXISTS "${NEP_LIBRARY}")
    message(FATAL_ERROR "NEP library not found under ${NEP_DIR}")
  endif()
endif()

add_executable(nep_benchmark nep_benchmark.cpp)
target_include_directories(nep_benchmark PRIVATE ${NEP_INCLUDE_DIR})
target_link_libraries(nep_benchmark PRIVATE ${NEP_LIBRARY})

if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(nep_benchmark PRIVATE OpenMP::OpenMP_CXX)
endif()
