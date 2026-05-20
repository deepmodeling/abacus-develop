# =============================================================================
# Setup Testing Environment (GTest, CTest, AddTest function)
# ==============================================================================

# include_guard(GLOBAL)

# --- Helper Macro: Ensure a minimum C++ standard version ---
macro(set_if_higher VARIABLE VALUE)
  if(${VARIABLE} LESS ${VALUE})
    set(${VARIABLE} ${VALUE})
  endif()
endmacro()

# Add performance test in abacus
if(ENABLE_GOOGLEBENCH)
  set(BUILD_TESTING ON)
  find_package(benchmark HINTS ${BENCHMARK_DIR})
  if(NOT ${benchmark_FOUND})
    set(BENCHMARK_USE_BUNDLED_GTEST OFF)
    include(FetchContent)
    FetchContent_Declare(
      benchmark
      GIT_REPOSITORY https://github.com/google/benchmark.git
      GIT_TAG "origin/main"
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE)
    set(BENCHMARK_ENABLE_TESTING OFF)
    FetchContent_MakeAvailable(benchmark)
  endif()
endif()

 function(AddTest) # function for UT
    cmake_parse_arguments(UT "DYN" "TARGET"
                          "LIBS;DYN_LIBS;STATIC_LIBS;SOURCES;DEPENDS" ${ARGN})
    add_executable(${UT_TARGET} ${UT_SOURCES})

    if(ENABLE_COVERAGE)
      add_coverage(${UT_TARGET})
    endif()

    # dependencies & link library
    target_link_libraries(${UT_TARGET} ${UT_LIBS} Threads::Threads
                          GTest::gtest_main GTest::gmock_main)
    if(ENABLE_GOOGLEBENCH)
      target_link_libraries(
        ${UT_TARGET} benchmark::benchmark)
    endif()

    if(USE_OPENMP)
      target_link_libraries(${UT_TARGET} OpenMP::OpenMP_CXX)
    endif()

    # Link to build info if needed
    if("${UT_SOURCES}" MATCHES "parse_args.cpp")
        target_include_directories(${UT_TARGET} PUBLIC ${CMAKE_BINARY_DIR}/source/source_io)
    endif()
        
    install(TARGETS ${UT_TARGET} DESTINATION ${CMAKE_BINARY_DIR}/tests)
    add_test(
      NAME ${UT_TARGET}
      COMMAND ${UT_TARGET}
      WORKING_DIRECTORY $<TARGET_FILE_DIR:${UT_TARGET}>)
  endfunction(AddTest)

if(BUILD_TESTING)
  set_if_higher(CMAKE_CXX_STANDARD 14) # Required in orbital
  include(CTest)
  enable_testing()
  find_package(GTest HINTS /usr/local/lib/ ${GTEST_DIR})
  if(NOT ${GTest_FOUND})
    include(FetchContent)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG "origin/main"
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE)
    FetchContent_MakeAvailable(googletest)
  endif()
  # TODO: Try the GoogleTest module.
  # https://cmake.org/cmake/help/latest/module/GoogleTest.html
  add_subdirectory(tests) # Contains integration tests
endif()
