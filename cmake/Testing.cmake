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

# --- Main function to configure everything related to testing ---
function(setup_testing)
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
    
    if(NOT BUILD_TESTING)
        message(STATUS "Testing is disabled via BUILD_TESTING=OFF.")
        return()
    endif()

    message(STATUS "Setting up testing environment...")

    # Set a minimum C++ standard for orbital
    set_if_higher(CMAKE_CXX_STANDARD 14)

    include(CTest)
    enable_testing()

    find_package(GTest QUIET HINTS /usr/local/lib/ ${GTEST_DIR})
    if(NOT GTest_FOUND)
        message(STATUS "GTest not found. Fetching from GitHub...")
        include(FetchContent)
        
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG "origin/main" # Consider pinning to a stable tag
            GIT_SHALLOW TRUE
            GIT_PROGRESS TRUE
        )
        # Prevent GTest from overriding our project's compiler flags
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()

    # Define the AddTest function for creating unit tests
    function(AddTest)
        cmake_parse_arguments(UT "DYN" "TARGET"
                              "LIBS;DYN_LIBS;STATIC_LIBS;SOURCES;DEPENDS" ${ARGN})
        add_executable(${UT_TARGET} ${UT_SOURCES})

        if(ENABLE_COVERAGE)
            add_coverage(${UT_TARGET})
        endif()

        target_link_libraries(${UT_TARGET} ${UT_LIBS} Threads::Threads
                              GTest::gtest_main GTest::gmock_main)
                              
        if(ENABLE_GOOGLEBENCH)
            target_link_libraries(${UT_TARGET} benchmark::benchmark)
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
    endfunction()

    # Add the test subdirectory
    if(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
        add_subdirectory(tests)
    elseif(EXISTS ${CMAKE_SOURCE_DIR}/test/CMakeLists.txt)
        add_subdirectory(test)
    else()
        message(WARNING "BUILD_TESTING is ON, but no 'tests/' or 'test/' directory found.")
    endif()

    message(STATUS "Testing environment configured successfully.")

endfunction()
