cmake_minimum_required(VERSION 3.16)
if(POLICY CMP0135) # https://cmake.org/cmake/help/git-stage/policy/CMP0135.html
  cmake_policy(SET CMP0135 NEW)
  # Otherwise this policy generates a warning on CMake 3.24
endif()
if(POLICY CMP0144) # https://cmake.org/cmake/help/git-stage/policy/CMP0144.html
  cmake_policy(SET CMP0144 NEW)
endif()

project(
  ABACUS
  DESCRIPTION "ABACUS is an electronic structure package based on DFT."
  HOMEPAGE_URL "https://github.com/deepmodeling/abacus-develop"
  LANGUAGES CXX)

option(ENABLE_MPI "Enable MPI" ON)
option(USE_OPENMP "Enable OpenMP" ON)
option(USE_CUDA "Enable CUDA" OFF)
option(USE_CUDA_MPI "Enable CUDA-aware MPI" OFF)
option(USE_CUDA_ON_DCU "Enable CUDA on DCU" OFF)
option(USE_ROCM "Enable ROCm" OFF)
option(USE_DSP "Enable DSP" OFF)
option(USE_KML "Enable Kunpeng Math Library" OFF)
option(USE_SW "Enable SW Architecture" OFF)

option(USE_ABACUS_LIBM "Build libmath from source to speed up" OFF)
option(ENABLE_LIBXC "Enable using the LibXC package" OFF)
option(ENABLE_FLOAT_FFTW "Enable using single-precision FFTW library." OFF)

# option(ENABLE_DEEPKS "Enable the DeePKS algorithm" OFF)
# option(ENABLE_MLKEDF "Enable the Machine-Learning-based KEDF for OFDFT" OFF)

option(ENABLE_MLALGO "Enable the machine learning algorithms" OFF)

option(ENABLE_LCAO "Enable LCAO algorithm" ON)
option(USE_ELPA "Enable ELPA for LCAO" ON)
option(ENABLE_LIBRI "Enable LibRI for hybrid functional" OFF)
option(ENABLE_LIBCOMM "Enable LibComm" OFF)
option(ENABLE_PEXSI "Enable PEXSI for LCAO" OFF)
option(ENABLE_DFTD4 "Enable DFT-D4 dispersion correction" OFF)

option(BUILD_TESTING "Build unittests" OFF)
option(DEBUG_INFO "Print message to debug" OFF)
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(INFO "Enable gathering math library information" OFF)
option(ENABLE_COVERAGE "Enable coverage build" OFF)
option(GIT_SUBMODULE "Check submodules during build" ON)

# Do not enable it if generated code will run on different CPUs
option(ENABLE_NATIVE_OPTIMIZATION
       "Enable compilation optimization for the native machine's CPU type" OFF)

option(COMMIT_INFO "Print commit information in log" ON)
option(ENABLE_FFT_TWO_CENTER "Enable FFT-based two-center integral method" ON)
option(ENABLE_GOOGLEBENCH "Enable GOOGLE-benchmark usage" OFF)
option(ENABLE_RAPIDJSON "Enable rapid-json usage" OFF)
option(ENABLE_CNPY "Enable cnpy usage" OFF)
option(ENABLE_CUSOLVERMP "Enable cusolvermp" OFF)
option(ENABLE_NCCL_PARALLEL_DEVICE "Enable NCCL-backed collectives in parallel_device" OFF)

# CTest defines BUILD_TESTING when it is first included. Include it only after
# ABACUS has declared its OFF-by-default option above.
include(CTest)

if(NOT DEFINED NVHPC_ROOT_DIR AND DEFINED ENV{NVHPC_ROOT})
  set(NVHPC_ROOT_DIR
      "$ENV{NVHPC_ROOT}"
      CACHE PATH "Path to NVIDIA HPC SDK root directory.")
endif()

# Feature definitions are collected while options and dependencies are resolved
# below. They are applied to targets in source/CMakeLists.txt.
set_property(GLOBAL PROPERTY ABACUS_FEATURE_DEFINITIONS "")

function(abacus_normalize_definitions out_var)
  set(_defs)
  foreach(_def IN LISTS ARGN)
    if(_def MATCHES "^-D(.+)")
      list(APPEND _defs "${CMAKE_MATCH_1}")
    else()
      list(APPEND _defs "${_def}")
    endif()
  endforeach()
  set(${out_var} ${_defs} PARENT_SCOPE)
endfunction()

function(abacus_add_feature_definitions)
  abacus_normalize_definitions(_defs ${ARGN})
  set_property(GLOBAL APPEND PROPERTY ABACUS_FEATURE_DEFINITIONS ${_defs})
endfunction()

# enable json support
if(ENABLE_RAPIDJSON)
  find_package(RapidJSON CONFIG REQUIRED)
  abacus_add_feature_definitions(__RAPIDJSON)
endif()

# get commit info
if(COMMIT_INFO)
  find_package(Git)
  if(NOT Git_FOUND)
    message(
      WARNING
        "Git is not found, and abacus will not output the git commit information in log. \n\
You can install Git first and reinstall abacus.")
  else()
    message(STATUS "Found git: attempting to get commit info...")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%h
      OUTPUT_VARIABLE GIT_COMMIT_HASH
      OUTPUT_STRIP_TRAILING_WHITESPACE
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      RESULT_VARIABLE GIT_COMMIT_HASH_RESULT)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} log -1 --format=%cd
      OUTPUT_VARIABLE GIT_COMMIT_DATE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      RESULT_VARIABLE GIT_COMMIT_DATE_RESULT)
    if(GIT_COMMIT_HASH_RESULT EQUAL 0 AND GIT_COMMIT_DATE_RESULT EQUAL 0)
      abacus_add_feature_definitions(COMMIT_INFO)
      file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/commit.h"
           "#define COMMIT \"${GIT_COMMIT_HASH} (${GIT_COMMIT_DATE})\"\n")
      set(ABACUS_COMMIT_INFO_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}")
      message(STATUS "Current commit hash: ${GIT_COMMIT_HASH}")
      message(STATUS "Last commit date: ${GIT_COMMIT_DATE}")
    else()
      message(WARNING "Failed to get git commit info")
    endif()
  endif()
endif()

# Serial version of ABACUS will not use ELPA
if(NOT ENABLE_MPI)
  set(USE_ELPA OFF)
  set(ENABLE_MLALGO OFF)
endif()

# Different exe files of ABACUS
unset(ABACUS_BIN_NAME CACHE)

# Case : LCAO or PW
if(ENABLE_LCAO)
  # Case: CUDA is enabled (suffix with 'g' for GPU)
  if(USE_CUDA)
    if(ENABLE_MPI)
      if(NOT ENABLE_LIBRI AND NOT ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_basic_gpu)
      elseif(ENABLE_LIBRI AND NOT ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_std_gpu)
      elseif(NOT ENABLE_LIBRI AND ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_ml_gpu)
      elseif(ENABLE_LIBRI AND ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_max_gpu)
      endif()
    endif()
  # Case: CPU is enabled (suffix with 'p' for parallel) 
  else()
    if(ENABLE_MPI)
      if(NOT ENABLE_LIBRI AND NOT ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_basic_para)
      elseif(ENABLE_LIBRI AND NOT ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_std_para)
      elseif(NOT ENABLE_LIBRI AND ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_ml_para)
      elseif(ENABLE_LIBRI AND ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_max_para)
      endif()
    else()
      if(NOT ENABLE_LIBRI AND NOT ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_basic_ser)
      elseif(ENABLE_LIBRI AND NOT ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_std_ser)
      elseif(NOT ENABLE_LIBRI AND ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_ml_ser)
      elseif(ENABLE_LIBRI AND ENABLE_MLALGO)
        set(ABACUS_BIN_NAME abacus_max_ser)
      endif()
    endif()
  endif()
# Case : PW only
else()
  if(USE_CUDA)
    if(ENABLE_MPI)
      set(ABACUS_BIN_NAME abacus_pw_gpu)
    endif()
  else()
    if(ENABLE_MPI)
      set(ABACUS_BIN_NAME abacus_pw_para)
    else()
      set(ABACUS_BIN_NAME abacus_pw_ser)
    endif()
  endif()
endif()


# Use DSP hardware
if (USE_DSP)
  set(USE_ELPA OFF)
  set(ABACUS_BIN_NAME abacus_dsp)
endif()

if (USE_CUDA_ON_DCU)
  abacus_add_feature_definitions(__CUDA_ON_DCU)
endif()

if (USE_CUDA_MPI)
  abacus_add_feature_definitions(__CUDA_MPI)
endif()

list(APPEND CMAKE_MODULE_PATH
     "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
     "${PROJECT_SOURCE_DIR}/cmake/modules")

if(ENABLE_COVERAGE)
  find_package(codecov)
  if(NOT codecov_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      cmakecodecov
      GIT_REPOSITORY https://github.com/RWTH-HPC/CMake-codecov.git
      GIT_TAG master
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE
    )
    FetchContent_Populate(cmakecodecov)
    list(APPEND CMAKE_MODULE_PATH ${cmakecodecov_SOURCE_DIR}/cmake)
    set(CMAKE_MODULE_PATH "${cmakecodecov_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
    find_package(codecov REQUIRED)
  endif()
endif()

set(ABACUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(ENABLE_DFTD4)
  # DFTD4 requires enabling C and Fortran to work
  enable_language(C)
  enable_language(Fortran)
  # Avoid custom-lapack fallback when resolving DFT-D4 dependencies
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(dftd4 4.2.0 REQUIRED)
endif()

macro(set_if_higher VARIABLE VALUE)
  if(${VARIABLE} LESS ${VALUE})
    set(${VARIABLE} ${VALUE})
  endif()
endmacro()
if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
  message(WARNING "GCC4 is not fully supported.")
endif()
set(FETCHCONTENT_QUIET FALSE) # Notify user when cloning git repo

find_program(CCACHE ccache)
if(CCACHE)
  set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE} ${CMAKE_CXX_COMPILER_LAUNCHER})
  set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE} ${CMAKE_C_COMPILER_LAUNCHER})
  set(CMAKE_CUDA_COMPILER_LAUNCHER ${CCACHE} ${CMAKE_CUDA_COMPILER_LAUNCHER})
endif()

# Choose build type from: Debug Release RelWithDebInfo MinSizeRel Select
# 'Release' configuration for best performance; this will disable all
# assertions. Other default configurations are also available, see:
# https://cmake.org/cmake/help/latest/manual/cmake-buildsystem.7.html#default-and-custom-configurations
# For default flags, see:
# https://github.com/Kitware/CMake/blob/master/Modules/Compiler/GNU.cmake#L55

if(ENABLE_COVERAGE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()
if(ENABLE_ASAN)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
  add_compile_options(-O3 -g)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  # stick to strict floating point model on Intel Compiler
  add_compile_options(-fp-model=strict)
  set(USE_ABACUS_LIBM OFF) # Force turn off USE_ABACUS_LIBM on Intel Compiler
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wno-write-strings "
  )
endif()

if(USE_ABACUS_LIBM)
  abacus_add_feature_definitions(USE_ABACUS_LIBM)
endif()

if(ENABLE_NATIVE_OPTIMIZATION)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")
endif()

# Windows (native build, e.g. MinGW-w64 or MSVC) portability defines:
#   _USE_MATH_DEFINES     - expose M_PI and friends from <cmath>
#   NOMINMAX              - stop <windows.h> defining min()/max() macros
#   _CRT_SECURE_NO_WARNINGS - silence CRT "use _s function" deprecations
if(WIN32)
  abacus_add_feature_definitions(_USE_MATH_DEFINES NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()

if(ENABLE_LCAO)
  find_package(cereal CONFIG REQUIRED)
  abacus_add_feature_definitions(__LCAO)
  if(USE_ELPA)
    find_package(ELPA REQUIRED)
    abacus_add_feature_definitions(__ELPA)
  endif()

  if(ENABLE_FFT_TWO_CENTER)
    abacus_add_feature_definitions(USE_NEW_TWO_CENTER)
  endif()

  if(ENABLE_PEXSI)
    find_package(PEXSI REQUIRED)
    abacus_add_feature_definitions(__PEXSI)
    set(CMAKE_CXX_STANDARD 14)
  endif()
else()
  set(ENABLE_MLALGO OFF)
  set(ENABLE_LIBRI OFF)
endif()

if(DEBUG_INFO)
  abacus_add_feature_definitions(__DEBUG)
endif()

if(ENABLE_MPI)
  find_package(MPI COMPONENTS CXX REQUIRED)
  abacus_add_feature_definitions(__MPI)
endif()


if (USE_DSP)
  abacus_add_feature_definitions(__DSP)
endif()



if (USE_SW)
  abacus_add_feature_definitions(__SW)
  set(SW ON)
endif()

find_package(Threads REQUIRED)

if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  if(NOT DEFINED USE_CUDA)
    message(
      "CUDA components detected. \nWill build the CUDA version of ABACUS by default."
    )
    set(USE_CUDA ON)
  else()
    if(NOT USE_CUDA)
      message(
        STATUS
          "CUDA components detected, but USE_CUDA is set to OFF. NOT building CUDA version of ABACUS."
      )
    endif()
  endif()
else() # CUDA not found
  if(USE_CUDA)
    message(
      FATAL_ERROR
        "USE_CUDA is set but no CUDA components found.")
  endif()
endif()

if(USE_CUDA)
  cmake_minimum_required(VERSION 3.18) # required by `CUDA_ARCHITECTURES` below
  
  # Always find CUDAToolkit to get CUDAToolkit_VERSION
  find_package(CUDAToolkit REQUIRED)
  
  set_if_higher(CMAKE_CXX_STANDARD 14)
  if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL "13.0")
    message(STATUS "CUDA ${CUDAToolkit_VERSION} detected. Setting CMAKE_CUDA_STANDARD to 17.")
    set_if_higher(CMAKE_CXX_STANDARD 17)
  endif()
  set(CMAKE_CXX_EXTENSIONS ON)
  set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    # check
    # https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/Internal/CMakeCUDAArchitecturesAll.cmake
    # for available architectures in different CUDA versions
    
    # CUDA 13.0+ dropped support for architectures below 75
    if(CUDAToolkit_VERSION VERSION_LESS "13.0")
      set(CMAKE_CUDA_ARCHITECTURES
          60 # P100
          70 # V100
      )
    else()
      # Start with empty list; architectures 75+ will be added below
      set(CMAKE_CUDA_ARCHITECTURES)
    endif()
    
    # Add your CUDA arch here Check the Compute Capability version of your
    # GPU at: https://en.wikipedia.org/wiki/CUDA#GPUs_supported
    
    if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 10.0)
      list(APPEND CMAKE_CUDA_ARCHITECTURES 75) # T4
    endif()
    if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 11.0)
      list(APPEND CMAKE_CUDA_ARCHITECTURES 80) # A100
    endif()
    if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 11.1)
      list(APPEND CMAKE_CUDA_ARCHITECTURES 86)
    endif()
    if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 11.8)
      list(APPEND CMAKE_CUDA_ARCHITECTURES 89 90)
    endif()
  endif()
  enable_language(CUDA)
  if(USE_CUDA)
    abacus_add_feature_definitions(__CUDA)
    abacus_add_feature_definitions(__UT_USE_CUDA)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -g -G" CACHE STRING "CUDA flags for debug build" FORCE)
    endif()
    if (USE_OPENMP AND OpenMP_CXX_FOUND)
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=${OpenMP_CXX_FLAGS}" CACHE STRING "CUDA flags" FORCE)
    endif()
    if (ENABLE_NCCL_PARALLEL_DEVICE)
      if (NOT ENABLE_MPI)
        message(FATAL_ERROR
          "ENABLE_NCCL_PARALLEL_DEVICE requires ENABLE_MPI=ON.")
      endif()
      abacus_add_feature_definitions(__NCCL_PARALLEL_DEVICE)
      include(cmake/modules/SetupNccl.cmake)
      abacus_setup_nccl()
    endif()
    if (ENABLE_CUSOLVERMP)
      # Keep cuSOLVERMp discovery/linking logic in a dedicated module.
      include(cmake/modules/SetupCuSolverMp.cmake)
      abacus_setup_cusolvermp()
    endif()
    if (ENABLE_CUBLASMP)
      # Enforcement 1: cuBLASMp requires cuSOLVERMp to be enabled
      if (NOT ENABLE_CUSOLVERMP)
        message(FATAL_ERROR 
          "ENABLE_CUBLASMP is set to ON, but ENABLE_CUSOLVERMP is OFF."
          "In ABACUS, cuBLASMp support requires cuSOLVERMp to be enabled simultaneously."
          "Please set -DENABLE_CUSOLVERMP=ON.")
      endif()
      # Enforcement 2: cuBLASMp 0.8.0+ is incompatible with CAL backend
      # Note: _use_cal is defined inside abacus_setup_cusolvermp
      if (_use_cal)
        message(FATAL_ERROR 
          "cuBLASMp 0.8.0+ requires NCCL Symmetric Memory, but cuSOLVERMp is using CAL backend."
          "Please upgrade cuSOLVERMp to >= 0.7.0 to use NCCL for both.")
      endif()
      include(cmake/modules/SetupCuBlasMp.cmake)
      abacus_setup_cublasmp()
    endif()
  endif()
endif()

# Warning: CMake add support to HIP in version 3.21. This is rather a new
# version. Use cmake with AMD-ROCm:
# https://rocmdocs.amd.com/en/latest/Installation_Guide/Using-CMake-with-AMD-ROCm.html
if(USE_ROCM)
  cmake_minimum_required(VERSION 3.18)
  if(NOT DEFINED ROCM_PATH)
    set(ROCM_PATH
        "/opt/rocm"
        CACHE STRING "Default ROCM installation directory.")
  endif()
  if(NOT DEFINED HIP_PATH)
    if(NOT DEFINED ENV{HIP_PATH})
      set(HIP_PATH
          "${ROCM_PATH}/hip"
          CACHE PATH "Path to which HIP has been installed")
    else()
      set(HIP_PATH
          $ENV{HIP_PATH}
          CACHE PATH "Path to which HIP has been installed")
    endif()
  endif()
  set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
  set(HIP_HIPCC_FLAGS -fno-gpu-rdc; --std=c++14) # --amdgpu-target=gfx906

  # find_package(ROCM REQUIRED)
  find_package(HIP REQUIRED)
  find_package(hipfft REQUIRED)
  find_package(hipblas REQUIRED)
  find_package(hipsolver REQUIRED)

  if(HIP_FOUND)
    message(STATUS "Found HIP: " ${HIP_VERSION})
  else()
    message(
      FATAL_ERROR
        "Could not find HIP. Ensure that HIP is either installed in ${ROCM_PATH}/hip or the variable HIP_PATH is set to point to the right location."
    )
  endif()

  abacus_add_feature_definitions(__ROCM)
  abacus_add_feature_definitions(__UT_USE_ROCM)
  abacus_add_feature_definitions(__HIP_PLATFORM_HCC__)
endif()

if(ENABLE_ASAN)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(
      FATAL_ERROR
        "Address Sanitizer is not supported on Intel Classic compiler. Use Intel oneAPI compiler (icpx) instead."
    )
  endif()
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address)
endif()

if(DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
  set(MKLROOT "$ENV{MKLROOT}")
endif()
if(USE_KML)
  set(_kml_components BLAS LAPACK FFTW3)
  if(ENABLE_MPI)
    list(APPEND _kml_components ScaLAPACK)
  endif()
  if(ENABLE_FLOAT_FFTW)
    list(APPEND _kml_components FFTW3_FLOAT)
  endif()

  find_package(KML REQUIRED COMPONENTS ${_kml_components})
  abacus_add_feature_definitions(__KML)
elseif(MKLROOT)
  set(MKL_INTERFACE lp64)
  set(ENABLE_SCALAPACK ON)
  find_package(MKL REQUIRED)
  abacus_add_feature_definitions(__MKL)
elseif(NOT USE_SW)
  find_package(Lapack REQUIRED)
  # ScaLAPACK is a distributed-memory library and is only needed for the
  # MPI build. A serial build (e.g. the native Windows serial version)
  # must not require it.
  if(ENABLE_MPI)
    find_package(ScaLAPACK REQUIRED)
  endif()
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Intel|Clang")
    message(WARNING "Cannot determine the required Fortran runtime.")
  endif()
endif()

if(NOT USE_KML AND NOT MKLROOT AND NOT USE_SW)
  find_package(FFTW3 REQUIRED)
endif()

if(ENABLE_FLOAT_FFTW)
  abacus_add_feature_definitions(__ENABLE_FLOAT_FFTW)
endif()

if(ENABLE_MLALGO)
  find_path(libnpy_SOURCE_DIR npy.hpp HINTS ${libnpy_INCLUDE_DIR})
  if(NOT libnpy_SOURCE_DIR)
    include(FetchContent)
    FetchContent_Declare(
      libnpy
      GIT_REPOSITORY https://github.com/llohse/libnpy.git
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE)
    FetchContent_MakeAvailable(libnpy)
  else()
  endif()

  abacus_add_feature_definitions(__MLALGO)
endif()

# Torch uses outdated components to detect CUDA arch, causing failure on
# latest CUDA kits. Set CMake variable TORCH_CUDA_ARCH_LIST in the form of
# "major.minor" if required.
if(ENABLE_MLALGO OR DEFINED Torch_DIR)
  find_package(Torch REQUIRED)
  if(NOT Torch_VERSION VERSION_LESS "2.1.0")
    set_if_higher(CMAKE_CXX_STANDARD 17)
  elseif(NOT Torch_VERSION VERSION_LESS "1.5.0")
    set_if_higher(CMAKE_CXX_STANDARD 14)
  endif()
  set(ABACUS_TORCH_CXX_FLAGS "${TORCH_CXX_FLAGS}")
endif()

if (ENABLE_CNPY)
  find_path(cnpy_SOURCE_DIR
    cnpy.h
    HINTS ${libnpy_INCLUDE_DIR}
  )
  if(NOT cnpy_SOURCE_DIR)
    include(FetchContent)
    FetchContent_Declare(
      cnpy
      GIT_REPOSITORY https://github.com/rogersce/cnpy.git
      GIT_PROGRESS TRUE
    )
    FetchContent_MakeAvailable(cnpy)
  else()
  endif()

  # find ZLIB and link
  find_package(ZLIB REQUIRED)
  abacus_add_feature_definitions(__USECNPY)
endif()

function(git_submodule_update)
  if(GIT_SUBMOD_RESULT EQUAL "0")
    message(DEBUG "Submodule init'ed")
  else()
    find_package(Git REQUIRED)
    if(NOT EXISTS "${PROJECT_SOURCE_DIR}/.git")
      message(FATAL_ERROR "Not a git repository; cannot update submodules")
    endif()

    message(STATUS "Updating submodules")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
          "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}."
      )
    endif()
  endif()
endfunction()

if(DEFINED LIBRI_DIR)
  set(ENABLE_LIBRI ON)
endif()
if(ENABLE_LIBRI)
  set_if_higher(CMAKE_CXX_STANDARD 14)
  if(LIBRI_DIR)
  else()
    find_package(LibRI REQUIRED)
  endif()
  abacus_add_feature_definitions(__EXX EXX_DM=3 EXX_H_COMM=2 TEST_EXX_LCAO=0
                                 TEST_EXX_RADIAL=1)
endif()

if(ENABLE_LIBRI OR DEFINED LIBCOMM_DIR)
  set(ENABLE_LIBCOMM ON)
endif()
if(ENABLE_LIBCOMM)
  if(LIBCOMM_DIR)
  else()
    find_package(LibComm REQUIRED)
  endif()
endif()


if(ENABLE_LIBXC)
  find_package(Libxc CONFIG REQUIRED)
  if(Libxc_VERSION VERSION_LESS "5.1.7")
    message(FATAL_ERROR "Libxc >= 5.1.7 is required")
  endif()
  abacus_add_feature_definitions(USE_LIBXC)
endif()

if(DEFINED DeePMD_DIR)
  abacus_add_feature_definitions(__DPMD HIGH_PREC)
  add_compile_options(-Wl,--no-as-needed)
  find_package(DeePMD REQUIRED)
  if(DeePMDC_FOUND)
    abacus_add_feature_definitions(__DPMDC)
  endif()
endif()

if(DEFINED NEP_DIR)
  find_package(NEP REQUIRED)

  if(NEP_FOUND)
    abacus_add_feature_definitions(__NEP)
  endif()
endif()

if(DEFINED TensorFlow_DIR)
  find_package(TensorFlow REQUIRED)
endif()

abacus_add_feature_definitions(__FFTW3 __SELINV METIS)

if(INFO)
  message(STATUS "Will gather math lib info.")
  abacus_add_feature_definitions(GATHER_INFO)
  # modifications on blas_connector and lapack_connector
endif()

# Target construction, final link ordering, installation, and test targets are
# intentionally centralized in source/CMakeLists.txt.
add_subdirectory(source)
