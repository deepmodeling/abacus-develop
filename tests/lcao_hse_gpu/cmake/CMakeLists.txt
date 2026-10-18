cmake_minimum_required(VERSION 3.18)
project(libri_consumer_check LANGUAGES CXX)

# Exercise the dependency target without ABACUS's global includes or link flags.
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../../../cmake/modules")
option(ENABLE_OPENMP "Enable OpenMP" ON)
option(ENABLE_LIBRI_GPU "Check the LibRI CUDA dependency interface" OFF)
option(ENABLE_LIBRI_GPU_TESTS "Run numerical LibRI tests requiring a GPU" OFF)
if(ENABLE_LIBRI_GPU_TESTS AND NOT ENABLE_LIBRI_GPU)
  message(FATAL_ERROR "ENABLE_LIBRI_GPU_TESTS requires ENABLE_LIBRI_GPU=ON")
endif()
find_package(MPI REQUIRED COMPONENTS CXX)
if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()
if(MKLROOT OR MKL_ROOT OR DEFINED ENV{MKLROOT})
  find_package(MKL REQUIRED)
else()
  find_package(Lapack REQUIRED)
endif()
if(ENABLE_LIBRI_GPU)
  find_package(CUDAToolkit REQUIRED)
endif()
find_package(LibRI REQUIRED)

add_executable(libri_consumer test_libri_consumer.cpp)
target_link_libraries(libri_consumer PRIVATE LibRI::LibRI)
target_compile_definitions(libri_consumer PRIVATE
  $<$<BOOL:${ENABLE_LIBRI_GPU}>:EXPECT_LIBRI_GPU>)
enable_testing()
add_test(NAME libri_consumer COMMAND libri_consumer)

if(ENABLE_LIBRI_GPU_TESTS)
  add_executable(libri_gpu_accumulation test_libri_gpu_accumulation.cpp)
  target_link_libraries(libri_gpu_accumulation PRIVATE LibRI::LibRI)
  add_test(NAME libri_gpu_accumulation
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:libri_gpu_accumulation> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(libri_gpu_accumulation PROPERTIES
    LABELS gpu
    TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
endif()
