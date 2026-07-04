# ==============================================================================
# External dependency interfaces
# ==============================================================================
#
# This file is the only place where ABACUS target linkage is assembled. The
# top-level CMakeLists.txt resolves options and packages but deliberately does
# not construct or link ABACUS targets.
#
# Link order is part of the build contract. Static linkers resolve archives
# from left to right, so every consumer must precede the provider of its
# unresolved symbols. Keep the final link closure ordered as follows:
#
#   1. ABACUS internal targets, from high-level consumers to lower-level
#      providers;
#   2. optional external feature libraries used by those internal targets;
#   3. numerical backends and their MPI, OpenMP, compiler-runtime, and system
#      dependencies.
#
# Do not sort or globally deduplicate this closure. A package target must carry
# its own transitive requirements in INTERFACE_LINK_LIBRARIES. Cyclic static
# archives must be handled locally by the corresponding package target or
# adapter rather than by wrapping the entire ABACUS link line in a linker group.

# This target intentionally carries compile usage requirements only. It must
# not link an external library: internal object libraries and unit tests consume
# it for headers, flags, and definitions, while the final executable owns the
# complete link closure below.
add_library(abacus_compile_requirements INTERFACE)
add_library(abacus::compile_requirements ALIAS abacus_compile_requirements)

target_include_directories(abacus_compile_requirements INTERFACE
  ${ABACUS_SOURCE_DIR}
  ${ABACUS_SOURCE_DIR}/source_base/module_container)

if(ABACUS_COMMIT_INFO_INCLUDE_DIR)
  target_include_directories(abacus_compile_requirements INTERFACE
    ${ABACUS_COMMIT_INFO_INCLUDE_DIR})
endif()

function(abacus_add_target_compile_requirements target)
  if(NOT TARGET "${target}")
    return()
  endif()

  get_target_property(_include_dirs "${target}" INTERFACE_INCLUDE_DIRECTORIES)
  if(_include_dirs)
    target_include_directories(abacus_compile_requirements INTERFACE
      ${_include_dirs})
  endif()

  get_target_property(_system_include_dirs "${target}" INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
  if(_system_include_dirs)
    target_include_directories(abacus_compile_requirements SYSTEM INTERFACE
      ${_system_include_dirs})
  endif()

  get_target_property(_compile_definitions "${target}" INTERFACE_COMPILE_DEFINITIONS)
  if(_compile_definitions)
    target_compile_definitions(abacus_compile_requirements INTERFACE
      ${_compile_definitions})
  endif()

  get_target_property(_compile_options "${target}" INTERFACE_COMPILE_OPTIONS)
  if(_compile_options)
    target_compile_options(abacus_compile_requirements INTERFACE
      ${_compile_options})
  endif()

  get_target_property(_compile_features "${target}" INTERFACE_COMPILE_FEATURES)
  if(_compile_features)
    target_compile_features(abacus_compile_requirements INTERFACE
      ${_compile_features})
  endif()
endfunction()

add_library(abacus_linalg_libs INTERFACE)
add_library(abacus::linalg_libs ALIAS abacus_linalg_libs)

add_library(abacus_feature_libs INTERFACE)
add_library(abacus::feature_libs ALIAS abacus_feature_libs)

add_library(abacus_link_libs INTERFACE)
add_library(abacus::link_libs ALIAS abacus_link_libs)

target_link_libraries(abacus_link_libs INTERFACE
  abacus::feature_libs
  abacus::linalg_libs)

# ------------------------------------------------------------------------------
# Ordered numerical backends and runtime closure
# ------------------------------------------------------------------------------

set(_abacus_linalg_libs)
set(_abacus_linalg_include_dirs)

if(USE_KML)
  if(ENABLE_MPI)
    list(APPEND _abacus_linalg_libs KML::ScaLAPACK)
  else()
    list(APPEND _abacus_linalg_libs KML::LAPACK)
  endif()
  list(APPEND _abacus_linalg_libs KML::FFTW3)
  list(APPEND _abacus_linalg_include_dirs ${KML_INCLUDE_DIRS})
  if(ENABLE_FLOAT_FFTW)
    list(APPEND _abacus_linalg_libs KML::FFTW3_FLOAT)
  endif()
elseif(MKLROOT)
  list(APPEND _abacus_linalg_libs MKL::MKL)
  list(APPEND _abacus_linalg_include_dirs ${MKL_INCLUDE} ${MKL_INCLUDE}/fftw)
  if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    list(APPEND _abacus_linalg_libs ifcore)
  endif()
elseif(USE_SW)
  list(APPEND _abacus_linalg_include_dirs
    ${SW_MATH}/include
    ${SW_FFT}/include)
  list(APPEND _abacus_linalg_libs
    ${SW_FFT}/lib/libfftw3.a
    ${SW_MATH}/libswfft.a
    ${SW_MATH}/libswscalapack.a
    ${SW_MATH}/libswlapack.a
    ${SW_MATH}/libswblas.a
    gfortran)
else()
  if(ENABLE_MPI)
    list(APPEND _abacus_linalg_libs ScaLAPACK::ScaLAPACK)
  endif()
  if(USE_OPENMP)
    list(APPEND _abacus_linalg_libs FFTW3::FFTW3_OMP)
  endif()
  list(APPEND _abacus_linalg_libs
    FFTW3::FFTW3
    LAPACK::LAPACK
    BLAS::BLAS)
  list(APPEND _abacus_linalg_include_dirs ${FFTW3_INCLUDE_DIRS})
  if(ENABLE_FLOAT_FFTW)
    list(APPEND _abacus_linalg_libs FFTW3::FFTW3_FLOAT)
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    list(APPEND _abacus_linalg_libs gfortran)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    list(APPEND _abacus_linalg_libs ifcore)
  endif()
endif()

if(ENABLE_MPI)
  list(APPEND _abacus_linalg_libs MPI::MPI_CXX)
endif()
if(USE_OPENMP)
  list(APPEND _abacus_linalg_libs OpenMP::OpenMP_CXX)
endif()
list(APPEND _abacus_linalg_libs Threads::Threads)

# libm exists on Linux and MinGW-w64 but not in the MSVC CRT.
if(NOT MSVC)
  list(APPEND _abacus_linalg_libs m)
endif()

target_link_libraries(abacus_linalg_libs INTERFACE ${_abacus_linalg_libs})
target_include_directories(abacus_linalg_libs INTERFACE ${_abacus_linalg_include_dirs})

target_include_directories(abacus_compile_requirements INTERFACE
  ${_abacus_linalg_include_dirs})

foreach(_abacus_linalg_target IN ITEMS
    Threads::Threads
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    FFTW3::FFTW3
    FFTW3::FFTW3_OMP
    FFTW3::FFTW3_FLOAT
    KML::LAPACK
    KML::ScaLAPACK
    KML::FFTW3
    KML::FFTW3_FLOAT)
  abacus_add_target_compile_requirements(${_abacus_linalg_target})
endforeach()

# ------------------------------------------------------------------------------
# Optional external feature libraries
# ------------------------------------------------------------------------------

set(_abacus_feature_libs)
set(_abacus_feature_include_dirs)
set(_abacus_feature_compile_options)

if(ENABLE_RAPIDJSON)
  list(APPEND _abacus_feature_libs RapidJSON)
endif()

if(ENABLE_LCAO)
  list(APPEND _abacus_feature_libs cereal::cereal)

  if(USE_ELPA)
    list(APPEND _abacus_feature_libs ELPA::ELPA)
    list(APPEND _abacus_feature_include_dirs ${ELPA_INCLUDE_DIR})
  endif()

  if(ENABLE_PEXSI)
    # Temporary adapter for the legacy FindPEXSI.cmake result. Replace this with
    # PEXSI::PEXSI when config-package discovery is adopted.
    list(APPEND _abacus_feature_libs
      ${PEXSI_LIBRARY}
      ${SuperLU_DIST_LIBRARY}
      ${ParMETIS_LIBRARY}
      ${METIS_LIBRARY})
    list(APPEND _abacus_feature_include_dirs
      ${PEXSI_INCLUDE_DIR}
      ${ParMETIS_INCLUDE_DIR})
  endif()
endif()

if(ENABLE_MLALGO)
  if(libnpy_INCLUDE_DIR)
    list(APPEND _abacus_feature_include_dirs ${libnpy_INCLUDE_DIR})
  endif()
  if(libnpy_SOURCE_DIR)
    list(APPEND _abacus_feature_include_dirs ${libnpy_SOURCE_DIR}/include)
  endif()
endif()

if(ENABLE_MLALGO OR DEFINED Torch_DIR)
  list(APPEND _abacus_feature_include_dirs ${TORCH_INCLUDE_DIRS})
  list(APPEND _abacus_feature_compile_options ${ABACUS_TORCH_CXX_FLAGS})
endif()

if(ENABLE_CNPY)
  list(APPEND _abacus_feature_libs cnpy ZLIB::ZLIB)
  if(cnpy_INCLUDE_DIR)
    list(APPEND _abacus_feature_include_dirs ${cnpy_INCLUDE_DIR})
  endif()
  if(cnpy_SOURCE_DIR)
    list(APPEND _abacus_feature_include_dirs ${cnpy_SOURCE_DIR})
  endif()
endif()

if(ENABLE_LIBRI)
  list(APPEND _abacus_feature_include_dirs ${LIBRI_DIR}/include)
endif()

if(ENABLE_LIBCOMM)
  list(APPEND _abacus_feature_include_dirs ${LIBCOMM_DIR}/include)
endif()

if(ENABLE_LIBXC)
  list(APPEND _abacus_feature_libs Libxc::xc)
endif()

if(DEFINED DeePMD_DIR)
  if(DeePMDC_FOUND)
    list(APPEND _abacus_feature_libs DeePMD::deepmd_c)
  else()
    list(APPEND _abacus_feature_libs DeePMD::deepmd_cc)
  endif()
  list(APPEND _abacus_feature_include_dirs ${DeePMD_DIR}/include)
endif()

if(DEFINED NEP_DIR AND NEP_FOUND)
  list(APPEND _abacus_feature_libs NEP::nep)
endif()

if(DEFINED TensorFlow_DIR AND TensorFlow_FOUND)
  list(APPEND _abacus_feature_libs TensorFlow::tensorflow_cc)
  list(APPEND _abacus_feature_include_dirs ${TensorFlow_DIR}/include)
endif()

if(USE_DSP)
  list(APPEND _abacus_feature_libs
    ${OMPI_LIBRARY1}
    ${MT_HOST_DIR}/hthreads/lib/libhthread_device.a
    ${MT_HOST_DIR}/hthreads/lib/libhthread_host.a)
  list(APPEND _abacus_feature_include_dirs
    ${MTBLAS_FFT_DIR}/libmtblas/include
    ${MT_HOST_DIR}/include)
endif()

if(USE_CUDA)
  if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 12.9)
    list(APPEND _abacus_feature_libs cudart)
  else()
    list(APPEND _abacus_feature_libs cudart nvToolsExt)
  endif()
  list(APPEND _abacus_feature_include_dirs ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 13.0
     AND EXISTS "${CUDAToolkit_ROOT}/include/cccl")
    list(APPEND _abacus_feature_include_dirs "${CUDAToolkit_ROOT}/include/cccl")
  endif()

  if(ENABLE_NCCL_PARALLEL_DEVICE)
    list(APPEND _abacus_feature_libs NCCL::NCCL)
  endif()
  if(ENABLE_CUSOLVERMP)
    if(_use_cal)
      list(APPEND _abacus_feature_libs CAL::CAL cusolverMp::cusolverMp)
    else()
      list(APPEND _abacus_feature_libs NCCL::NCCL cusolverMp::cusolverMp)
    endif()
  endif()
  if(ENABLE_CUBLASMP)
    list(APPEND _abacus_feature_libs cublasMp::cublasMp)
  endif()
endif()

if(USE_ROCM)
  list(APPEND _abacus_feature_libs
    hip::host
    hip::device
    hip::hipfft
    roc::hipblas
    roc::hipsolver)
  list(APPEND _abacus_feature_include_dirs ${ROCM_PATH}/include)
endif()

if(ENABLE_ASAN)
  list(APPEND _abacus_feature_libs -fsanitize=address)
endif()

target_link_libraries(abacus_feature_libs INTERFACE ${_abacus_feature_libs})

target_include_directories(abacus_compile_requirements INTERFACE
  ${_abacus_feature_include_dirs})
target_compile_options(abacus_compile_requirements INTERFACE
  ${_abacus_feature_compile_options})

foreach(_abacus_feature_target IN ITEMS
    RapidJSON
    cereal::cereal
    ELPA::ELPA
    Libxc::xc
    dftd4::dftd4
    DeePMD::deepmd_c
    DeePMD::deepmd_cc
    NEP::nep
    TensorFlow::tensorflow_cc
    ZLIB::ZLIB
    NCCL::NCCL
    CAL::CAL
    cusolverMp::cusolverMp
    cublasMp::cublasMp
    hip::host
    hip::device
    hip::hipfft
    roc::hipblas
    roc::hipsolver)
  abacus_add_target_compile_requirements(${_abacus_feature_target})
endforeach()

# ------------------------------------------------------------------------------
# Per-target feature definitions and common compile usage requirements
# ------------------------------------------------------------------------------

define_property(
  DIRECTORY
  PROPERTY ABACUS_DISABLED_FEATURE_DEFINITIONS
  INHERITED
  BRIEF_DOCS "ABACUS feature definitions disabled for targets in this directory"
  FULL_DOCS "Feature definitions disabled for targets created in this directory.")

define_property(
  DIRECTORY
  PROPERTY ABACUS_LOCAL_FEATURE_DEFINITIONS
  INHERITED
  BRIEF_DOCS "Additional ABACUS feature definitions for targets in this directory"
  FULL_DOCS "Additional feature definitions for targets created in this directory.")

function(abacus_disable_feature_definitions)
  abacus_normalize_definitions(_defs ${ARGN})
  set_property(DIRECTORY APPEND PROPERTY ABACUS_DISABLED_FEATURE_DEFINITIONS ${_defs})
endfunction()

function(abacus_add_local_feature_definitions)
  abacus_normalize_definitions(_defs ${ARGN})
  set_property(DIRECTORY APPEND PROPERTY ABACUS_LOCAL_FEATURE_DEFINITIONS ${_defs})
endfunction()

function(abacus_apply_build_options target)
  if(NOT TARGET "${target}")
    return()
  endif()

  get_target_property(_type "${target}" TYPE)
  if(_type STREQUAL "INTERFACE_LIBRARY" OR _type STREQUAL "UTILITY")
    return()
  endif()

  get_target_property(_imported "${target}" IMPORTED)
  if(_imported)
    return()
  endif()

  get_target_property(_source_dir "${target}" SOURCE_DIR)
  get_property(_defs GLOBAL PROPERTY ABACUS_FEATURE_DEFINITIONS)
  get_property(_disabled DIRECTORY "${_source_dir}" PROPERTY ABACUS_DISABLED_FEATURE_DEFINITIONS)
  get_property(_local DIRECTORY "${_source_dir}" PROPERTY ABACUS_LOCAL_FEATURE_DEFINITIONS)

  if(_disabled)
    list(REMOVE_ITEM _defs ${_disabled})
  endif()
  if(_local)
    list(APPEND _defs ${_local})
  endif()
  if(_defs)
    list(REMOVE_DUPLICATES _defs)
    target_compile_definitions("${target}" PRIVATE ${_defs})
  endif()

  target_link_libraries("${target}" PRIVATE abacus::compile_requirements)
endfunction()

function(abacus_apply_build_options_to_dir dir)
  # FetchContent adds third-party projects as subdirectories of its caller.
  # Keep the recursive application strictly within ABACUS's own source tree so
  # exported third-party targets never acquire ABACUS compile requirements.
  set(_abacus_source_root "${PROJECT_SOURCE_DIR}/source/")
  string(FIND "${dir}/" "${_abacus_source_root}" _abacus_source_pos)
  if(NOT _abacus_source_pos EQUAL 0)
    return()
  endif()

  get_property(_targets DIRECTORY "${dir}" PROPERTY BUILDSYSTEM_TARGETS)
  foreach(_target IN LISTS _targets)
    abacus_apply_build_options("${_target}")
  endforeach()

  get_property(_subdirs DIRECTORY "${dir}" PROPERTY SUBDIRECTORIES)
  foreach(_subdir IN LISTS _subdirs)
    abacus_apply_build_options_to_dir("${_subdir}")
  endforeach()
endfunction()

# ==============================================================================
# Executable, source tree, and tests
# ==============================================================================

set(ABACUS_TEST_DIR "${PROJECT_SOURCE_DIR}/tests")
include(${PROJECT_SOURCE_DIR}/cmake/Testing.cmake)

add_executable(${ABACUS_BIN_NAME} source_main/main.cpp)
set(ABACUS_BIN_PATH ${CMAKE_CURRENT_BINARY_DIR}/${ABACUS_BIN_NAME})

if(USE_CUDA)
  set_property(TARGET ${ABACUS_BIN_NAME}
    PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
  target_compile_definitions(${ABACUS_BIN_NAME} PRIVATE __USE_NVTX)
endif()

if(ENABLE_COVERAGE)
  add_coverage(${ABACUS_BIN_NAME})
endif()

include(${PROJECT_SOURCE_DIR}/cmake/BuildInfo.cmake)
setup_build_info()

add_subdirectory(source_base)
add_subdirectory(source_cell)
add_subdirectory(source_psi)
add_subdirectory(source_estate)
add_subdirectory(source_hamilt)
add_subdirectory(source_pw)
add_subdirectory(source_lcao)
add_subdirectory(source_hsolver)
add_subdirectory(source_basis/module_ao)
add_subdirectory(source_basis/module_nao)
add_subdirectory(source_md)
add_subdirectory(source_basis/module_pw)
add_subdirectory(source_esolver)
add_subdirectory(source_lcao/module_gint)
add_subdirectory(source_io)
add_subdirectory(source_relax)
add_subdirectory(source_lcao/module_ri)
add_subdirectory(source_io/module_parameter)
add_subdirectory(source_lcao/module_lr)
add_subdirectory(source_lcao/module_rdmft) # add by jghan

add_library(
  driver
  OBJECT
  source_main/driver.cpp
  source_main/driver_run.cpp)

list(APPEND device_srcs
  source_pw/module_pwdft/kernels/nonlocal_op.cpp
  source_pw/module_pwdft/kernels/veff_op.cpp
  source_pw/module_pwdft/kernels/ekinetic_op.cpp
  source_pw/module_pwdft/kernels/meta_op.cpp
  source_pw/module_stodft/kernels/hpsi_norm_op.cpp
  source_basis/module_pw/kernels/pw_op.cpp
  source_hsolver/kernels/hegvd_op.cpp
  source_hsolver/kernels/bpcg_kernel_op.cpp
  source_estate/kernels/elecstate_op.cpp

  source_base/module_device/device.cpp
  source_base/module_device/device_helpers.cpp
  source_base/module_device/output_device.cpp
  source_base/module_device/memory_op.cpp
  source_base/kernels/math_kernel_op.cpp
  source_base/kernels/math_kernel_op_vec.cpp

  source_pw/module_pwdft/kernels/force_op.cpp
  source_pw/module_pwdft/kernels/stress_op.cpp
  source_pw/module_pwdft/kernels/onsite_op.cpp
  source_pw/module_pwdft/kernels/wf_op.cpp
  source_pw/module_pwdft/kernels/vnl_op.cpp
  source_base/kernels/math_ylm_op.cpp
  source_hamilt/module_xc/kernels/xc_functional_op.cpp
  source_pw/module_pwdft/kernels/cal_density_real_op.cpp
  source_pw/module_pwdft/kernels/mul_potential_op.cpp
  source_pw/module_pwdft/kernels/vec_mul_vec_complex_op.cpp
  source_pw/module_pwdft/kernels/exx_cal_energy_op.cpp)

if(USE_CUDA)
  list(APPEND device_srcs
    source_pw/module_pwdft/kernels/cuda/nonlocal_op.cu
    source_pw/module_pwdft/kernels/cuda/veff_op.cu
    source_pw/module_pwdft/kernels/cuda/ekinetic_op.cu
    source_pw/module_pwdft/kernels/cuda/meta_op.cu
    source_pw/module_stodft/kernels/cuda/hpsi_norm_op.cu
    source_pw/module_pwdft/kernels/cuda/onsite_op.cu
    source_basis/module_pw/kernels/cuda/pw_op.cu
    source_hsolver/kernels/cuda/hegvd_op.cu
    source_hsolver/kernels/cuda/bpcg_kernel_op.cu
    source_estate/kernels/cuda/elecstate_op.cu
    source_base/module_device/cuda/memory_op.cu
    source_pw/module_pwdft/kernels/cuda/force_op.cu
    source_pw/module_pwdft/kernels/cuda/stress_op.cu
    source_pw/module_pwdft/kernels/cuda/wf_op.cu
    source_pw/module_pwdft/kernels/cuda/vnl_op.cu
    source_base/module_device/cuda_compat.cpp
    source_base/kernels/cuda/math_ylm_op.cu
    source_base/kernels/cuda/math_kernel_op.cu
    source_base/kernels/cuda/math_kernel_op_vec.cu
    source_hamilt/module_xc/kernels/cuda/xc_functional_op.cu
    source_pw/module_ofdft/kernels/cuda/kedf_wt_gpu.cu
    source_pw/module_pwdft/kernels/cuda/cal_density_real_op.cu
    source_pw/module_pwdft/kernels/cuda/mul_potential_op.cu
    source_pw/module_pwdft/kernels/cuda/vec_mul_vec_complex.cu
    source_pw/module_pwdft/kernels/cuda/exx_cal_energy_op.cu)
endif()

if(USE_ROCM)
  hip_add_library(device_rocm STATIC
    source_pw/module_pwdft/kernels/rocm/nonlocal_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/veff_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/ekinetic_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/meta_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/onsite_op.hip.cu
    source_pw/module_stodft/kernels/rocm/hpsi_norm_op.hip.cu
    source_basis/module_pw/kernels/rocm/pw_op.hip.cu
    source_hsolver/kernels/rocm/hegvd_op.hip.cu
    source_hsolver/kernels/rocm/bpcg_kernel_op.hip.cu
    source_estate/kernels/rocm/elecstate_op.hip.cu
    source_base/module_device/rocm/memory_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/force_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/stress_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/wf_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/vnl_op.hip.cu
    source_base/kernels/rocm/math_kernel_op.hip.cu
    source_base/kernels/rocm/math_kernel_op_vec.hip.cu
    source_base/kernels/rocm/math_ylm_op.hip.cu
    source_hamilt/module_xc/kernels/rocm/xc_functional_op.hip.cu)
endif()

if(USE_DSP)
  list(APPEND device_srcs
    source_base/kernels/dsp/dsp_connector.cpp)
endif()

add_library(device OBJECT ${device_srcs})

if(USE_CUDA)
  target_link_libraries(device PRIVATE cusolver cublas cufft)
elseif(USE_ROCM)
  target_link_libraries(device PRIVATE
    device_rocm
    hip::host
    hip::device
    hip::hipfft
    roc::hipblas
    roc::hipsolver)
endif()

# base uses symbols from device (memory_op and math_ylm_op).
target_link_libraries(base PUBLIC device)

if(ENABLE_COVERAGE)
  add_coverage(driver)
endif()

# ==============================================================================
# Final link closure and ordering
# ==============================================================================

target_link_libraries(
  ${ABACUS_BIN_NAME}
  PRIVATE
    # Internal ABACUS targets: consumers before providers.
    driver
    esolver
    hsolver
    hamilt_general
    elecstate
    module_pwdft
    module_ofdft
    module_stodft
    module_dfpt
    xc_
    vdw
    relax
    io_advanced
    io_basic
    io_input
    surchem
    neighbor_search
    neighbor
    md
    planewave
    symmetry
    cell
    parameter
    psi_overall_init
    psi_initializer
    psi
    dftu
    deltaspin
    container
    device
    base)

if(ENABLE_LCAO)
  target_link_libraries(
    ${ABACUS_BIN_NAME}
    PRIVATE
      hamilt_lcao
      tddft
      orb
      gint
      hcontainer
      numerical_atomic_orbitals
      lr
      rdmft)
  if(USE_ELPA)
    target_link_libraries(${ABACUS_BIN_NAME} PRIVATE genelpa)
  endif()
  if(ENABLE_PEXSI)
    target_link_libraries(${ABACUS_BIN_NAME} PRIVATE pexsi)
  endif()
  if(USE_CUDA)
    target_link_libraries(${ABACUS_BIN_NAME} PRIVATE diag_cusolver)
  endif()
endif()

if(ENABLE_MLALGO)
  target_link_libraries(${ABACUS_BIN_NAME} PRIVATE deepks hamilt_mlkedf)
endif()
if(ENABLE_LIBRI)
  target_link_libraries(${ABACUS_BIN_NAME} PRIVATE ri module_exx_symmetry)
endif()
if(ENABLE_RAPIDJSON)
  target_link_libraries(${ABACUS_BIN_NAME} PRIVATE json_output)
endif()

# External closure: feature consumers precede numerical providers and runtimes.
target_link_libraries(${ABACUS_BIN_NAME} PRIVATE abacus::link_libs)

# Torch includes some LAPACK routines, but with floating-point exceptions.
if(ENABLE_MLALGO OR DEFINED Torch_DIR)
  target_link_libraries(${ABACUS_BIN_NAME} PRIVATE ${TORCH_LIBRARIES})
endif()

# Apply feature definitions and common compile usage requirements only after all
# source targets have been created.
abacus_apply_build_options_to_dir("${CMAKE_CURRENT_SOURCE_DIR}")
abacus_apply_build_options(${ABACUS_BIN_NAME})

# Register integration tests only after the final executable and its path are
# available. Unit tests are added by source subdirectories through AddTest().
if(BUILD_TESTING)
  add_subdirectory("${ABACUS_TEST_DIR}" "${PROJECT_BINARY_DIR}/tests")
endif()

install(PROGRAMS ${ABACUS_BIN_PATH} TYPE BIN)

# Create a symbolic link 'abacus' pointing to the actual executable. Skipped on
# Windows because symlink creation needs elevated/developer-mode privileges and
# the executable carries an .exe suffix anyway.
if(NOT WIN32)
  install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${ABACUS_BIN_NAME} ${CMAKE_INSTALL_PREFIX}/bin/abacus WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)")
endif()

if(ENABLE_COVERAGE)
  coverage_evaluate()
endif()
