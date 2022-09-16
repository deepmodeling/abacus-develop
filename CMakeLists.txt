########################################
# CMake build system
# This file is part of ABACUS
cmake_minimum_required(VERSION 3.16)
########################################

project(ABACUS
    VERSION 2.3.0
    DESCRIPTION "ABACUS is an electronic structure package based on DFT."
    HOMEPAGE_URL "https://github.com/deepmodeling/abacus-develop"
    LANGUAGES CXX
)

option(ENABLE_DEEPKS "Enable DeePKS functionality" OFF)
option(ENABLE_LIBXC "Enable LibXC functionality" OFF)
option(USE_CUDA "Enable support to CUDA for PW." OFF)
option(USE_CUSOLVER_LCAO "Enable support to CUSOLVER for LCAO." OFF)
option(USE_ROCM "Enable support to ROCm." OFF)
option(USE_OPENMP " Enable OpenMP in abacus." ON)
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(BUILD_TESTING "Build ABACUS unit tests" OFF)
option(GENERATE_TEST_REPORTS "Enable test report generation" OFF)
option(INFO "Enable gathering of math library information" OFF)

set(ABACUS_BIN_NAME abacus)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)
set(ABACUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(ABACUS_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(ABACUS_BIN_PATH ${CMAKE_CURRENT_BINARY_DIR}/${ABACUS_BIN_NAME})
include_directories(${ABACUS_SOURCE_DIR})
add_executable(${ABACUS_BIN_NAME} source/main.cpp)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wno-write-strings)
set(FETCHCONTENT_QUIET FALSE) # Notify user when cloning git repo

find_program(CCACHE ccache)
if (CCACHE)
  set (CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE} ${CMAKE_CXX_COMPILER_LAUNCHER})
  set (CMAKE_C_COMPILER_LAUNCHER ${CCACHE} ${CMAKE_C_COMPILER_LAUNCHER})
endif()

# Select 'Release' configuration for best performance;
#   this will disable all assertions.
# Other default configurations are also available, see:
# https://cmake.org/cmake/help/latest/manual/cmake-buildsystem.7.html#default-and-custom-configurations
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  add_compile_options(-O2 -g) # default flag
endif()

find_package(Cereal REQUIRED)
include_directories(${CEREAL_INCLUDE_DIR})
add_compile_definitions(USE_CEREAL_SERIALIZATION)

find_package(ELPA REQUIRED)
include_directories(${ELPA_INCLUDE_DIR})
target_link_libraries(${ABACUS_BIN_NAME} ELPA::ELPA)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
target_link_libraries(${ABACUS_BIN_NAME} MPI::MPI_CXX)
add_compile_definitions(__MPI)
list(APPEND math_libs MPI::MPI_CXX)

find_package(Threads REQUIRED)
target_link_libraries(${ABACUS_BIN_NAME} Threads::Threads)

if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(${ABACUS_BIN_NAME} OpenMP::OpenMP_CXX)
  add_compile_options(${OpenMP_CXX_FLAGS})
  add_link_options(${OpenMP_CXX_LIBRARIES})
endif()

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  if(NOT DEFINED USE_CUDA OR NOT DEFINED USE_CUSOLVER_LCAO)
    if (NOT DEFINED USE_CUDA AND NOT DEFINED USE_CUSOLVER_LCAO)
      message("CUDA components detected. \nWill build the CUDA for PW version of ABACUS by default.")
      set(USE_CUDA ON)
      set(USE_CUSOLVER_LCAO OFF)
    elseif (NOT DEFINED USE_CUDA)
      set(USE_CUDA OFF)
    else()
      set(USE_CUSOLVER_LCAO OFF)
    endif()
  else()
    if(NOT USE_CUDA AND NOT USE_CUSOLVER_LCAO)
      message(STATUS "CUDA components detected, but both USE_CUDA and USE_CUSOLVER_LCAO set to OFF. NOT building CUDA version of ABACUS.")
    elseif (USE_CUDA AND USE_CUSOLVER_LCAO)
      message(FATAL_ERROR "USE_CUDA and USE_CUSOLVER_LCAO set, but now they not allowed to coexist.")
    endif()
  endif()
else() # CUDA not found
  if (USE_CUDA OR USE_CUSOLVER_LCAO)
    message(FATAL_ERROR "USE_CUDA or USE_CUSOLVER_LCAO set but no CUDA components found.")
    set(USE_CUDA OFF)
    set(USE_CUSOLVER_LCAO OFF)
  endif()
endif()

if(USE_CUDA OR USE_CUSOLVER_LCAO)
  cmake_minimum_required(VERSION 3.18) # required by `CUDA_ARCHITECTURES` below
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_EXTENSIONS ON)
  set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  enable_language(CUDA)
  target_link_libraries(${ABACUS_BIN_NAME}
    -lcufft
    -lcublas
    -lcudart
    -lcusolver
    -lnvToolsExt
  )
  set_property(TARGET ${ABACUS_BIN_NAME}
  PROPERTY CUDA_ARCHITECTURES
  60 # P100
  70 # V100
  75 # T4
  80 # A100
  )
  include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  if (USE_CUDA)
    add_compile_definitions(__CUDA)
  endif()
  if (USE_CUSOLVER_LCAO)
    add_compile_definitions(__CUSOLVER_LCAO)
  endif()
endif()


# Warning: CMake add support to HIP in version 3.21. This is rather a new version.
# Use cmake with AMD-ROCm: https://rocmdocs.amd.com/en/latest/Installation_Guide/Using-CMake-with-AMD-ROCm.html
if(USE_ROCM)
  SET(CMAKE_HIP_STANDARD ${CMAKE_CXX_STANDARD})
  list(APPEND CMAKE_PREFIX_PATH /opt/rocm/hip /opt/rocm)
  list(APPEND CMAKE_SYSTEM_PREFIX_PATH "/opt/rocm/hip")

  # find_package(ROCM REQUIRED)
  find_package(hip REQUIRED)
  find_package(hipfft REQUIRED)
  find_package(hipblas REQUIRED)
  include_directories(
    /opt/rocm/hip/include
    /opt/rocm/hipfft/include
    /opt/rocm/hipblas/include
  )
  target_link_libraries(${ABACUS_BIN_NAME}
    hip::device
    hip::host
    hip::hipfft
    roc::hipblas
  )
  add_compile_definitions(__ROCM)
endif()

if(ENABLE_ASAN)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(FATAL_ERROR "Address Sanitizer is not supported on Intel compiler.")
  endif()
  add_compile_options(
    -fsanitize=address
    -fno-omit-frame-pointer
  )
  add_link_options(
    -fsanitize=address
  )
  # `add_link_options` only affects executables added after.
  target_link_libraries(${ABACUS_BIN_NAME} -fsanitize=address)
endif()

if(DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
    set(MKLROOT "$ENV{MKLROOT}")
endif()
if(MKLROOT)
  find_package(IntelMKL REQUIRED)
  add_definitions(-D__MKL)
  include_directories(${MKL_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS}/fftw)

  # Since libtorch will find its own MKL, the fftw part conflicts with the original one.
  # When enable deepks, mkl will be linked within ${TORCH_LIBRARIES}.
  if(NOT ENABLE_DEEPKS)
    list(APPEND math_libs IntelMKL::MKL)
  endif()
else()
  find_package(FFTW3 REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(ScaLAPACK REQUIRED)
  include_directories(${FFTW3_INCLUDE_DIRS})
  list(APPEND math_libs
    FFTW3::FFTW3
    LAPACK::LAPACK
    ScaLAPACK::ScaLAPACK
  )
  if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    list(APPEND math_libs -lgfortran)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    list(APPEND math_libs -lifcore)
  else()
    message(WARNING "Cannot find the correct library for Fortran.")
  endif()
endif()

if(ENABLE_DEEPKS)
  set(CMAKE_CXX_STANDARD 14)
  find_package(Torch REQUIRED)
  include_directories(${TORCH_INCLUDE_DIRS})
  target_link_libraries(${ABACUS_BIN_NAME} deepks)
  list(APPEND math_libs ${TORCH_LIBRARIES})

  add_compile_options(${TORCH_CXX_FLAGS})

  find_path(libnpy_SOURCE_DIR
    npy.hpp
    HINTS ${libnpy_INCLUDE_DIR}
  )
  if(NOT libnpy_SOURCE_DIR)
    include(FetchContent)
    FetchContent_Declare(
      libnpy
      GIT_REPOSITORY https://github.com/llohse/libnpy.git
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE
    )
    FetchContent_MakeAvailable(libnpy)
  endif()
  include_directories(${libnpy_SOURCE_DIR}/include)
  add_compile_definitions(__DEEPKS)
endif()

list(APPEND math_libs m)
target_link_libraries(${ABACUS_BIN_NAME} ${math_libs})

if(DEFINED Libxc_DIR)
  set(ENABLE_LIBXC ON)
endif()
if(ENABLE_LIBXC)
  find_package(Libxc HINTS ${Libxc_DIR}/share/cmake/Libxc)
  if(${Libxc_FOUND})
    message("Using Libxc.")
    add_compile_definitions(USE_LIBXC)
    target_link_libraries(${ABACUS_BIN_NAME} Libxc::xc)
    include_directories(${Libxc_INCLUDE_DIRS})
  else()
    include(FetchContent)
    FetchContent_Declare(
      Libxc
      GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
      GIT_TAG "5.2.3"
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE
    )
    FetchContent_MakeAvailable(Libxc)
  endif()
endif()

add_compile_definitions(
  __EXX
  __FFTW3
  __FP
  __SELINV
  __LCAO
  METIS
  EXX_DM=3
  EXX_H_COMM=2
  TEST_EXX_LCAO=0
  TEST_EXX_RADIAL=1
)

if(INFO)
  message(STATUS "Will gather math lib info.")
  add_compile_definitions(GATHER_INFO)
  # modifications on blas_connector and lapack_connector
endif()

IF (BUILD_TESTING)
  set(CMAKE_CXX_STANDARD 14) # Required in orbital
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
      GIT_PROGRESS TRUE
    )
    FetchContent_MakeAvailable(googletest)
  endif()
  # TODO: Try the GoogleTest module.
  # https://cmake.org/cmake/help/latest/module/GoogleTest.html
  add_subdirectory(tests) # Contains integration tests

  function(AddTest) # function for UT
    cmake_parse_arguments(UT "DYN" "TARGET" "LIBS;DYN_LIBS;STATIC_LIBS;SOURCES;DEPENDS" ${ARGN})
    add_executable(${UT_TARGET} ${UT_SOURCES})
    #dependencies & link library
    target_link_libraries(${UT_TARGET} ${UT_LIBS}
      Threads::Threads GTest::gtest_main GTest::gmock_main)
    if(USE_OPENMP)
      target_link_libraries(${UT_TARGET} OpenMP::OpenMP_CXX)
    endif()
    install(TARGETS ${UT_TARGET} DESTINATION ${CMAKE_BINARY_DIR}/tests )
    add_test(NAME ${UT_TARGET}
      COMMAND ${UT_TARGET}
      WORKING_DIRECTORY $<TARGET_FILE_DIR:${UT_TARGET}>
    )
  endfunction(AddTest)
endif()

add_subdirectory(source)

target_link_libraries(${ABACUS_BIN_NAME}
    base
    cell
    symmetry
    md
    planewave
    surchem
    neighbor
    orb
    io
    relax
    lcao
    gint
    parallel
    mrrr
    pdiag
    pw
    ri
    driver
    xc_
    hsolver
    genelpa
    elecstate
    hamilt
    psi
    esolver
)

install(PROGRAMS ${ABACUS_BIN_PATH}
    TYPE BIN
    #DESTINATION ${CMAKE_INSTALL_BINDIR}
)
