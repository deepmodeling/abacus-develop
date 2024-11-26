cmake_minimum_required(VERSION 3.15...3.26)

# project settings
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

# find python and pybind11
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# set source path
set(ABACUS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../../source")
set(BASE_PATH "${ABACUS_SOURCE_DIR}/module_base")
set(NAO_PATH "${ABACUS_SOURCE_DIR}/module_basis/module_nao")
set(HSOLVER_PATH "${ABACUS_SOURCE_DIR}/module_hsolver")
set(PSI_PATH "${ABACUS_SOURCE_DIR}/module_psi")
set(ENABLE_LCAO ON)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake")

# add math_libs
if(DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
    set(MKLROOT "$ENV{MKLROOT}")
endif()
if(MKLROOT)
  set(MKL_INTERFACE lp64)
  set(ENABLE_MPI ON)
  if (ENABLE_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    #add_compile_definitions(__MPI)
  endif()

  set(USE_OPENMP ON)
  if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_link_options(${OpenMP_CXX_LIBRARIES})
  endif()
  find_package(MKL REQUIRED)
  add_definitions(-D__MKL)
  include_directories(${MKL_INCLUDE} ${MKL_INCLUDE}/fftw)

# Since libtorch will find its own MKL, the fftw part conflicts with the original one.
# When enable deepks, mkl will be linked within ${TORCH_LIBRARIES}.
  if(NOT ENABLE_DEEPKS)
    list(APPEND math_libs IntelMKL::MKL)
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    list(APPEND math_libs -lifcore)
  endif()
else()
    find_package(FFTW3 REQUIRED)
    add_compile_definitions(__FFTW3)
    find_package(LAPACK REQUIRED)
    include_directories(${FFTW3_INCLUDE_DIRS})
    message(STATUS "FFTW3_INCLUDE_DIRS: ${FFTW3_INCLUDE_DIRS}")
    list(APPEND math_libs FFTW3::FFTW3 LAPACK::LAPACK)

  if(ENABLE_LCAO)
    find_package(ScaLAPACK REQUIRED)
    list(APPEND math_libs ScaLAPACK::ScaLAPACK)
  endif()
endif()

# add include directories
include_directories(
    ${BASE_PATH} 
    ${ABACUS_SOURCE_DIR}
    ${ABACUS_SOURCE_DIR}/module_base/module_container
    )

# add basic libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# add base
set(BASE_BINARY_DIR "${PROJECT_SOURCE_DIR}/build/base")
add_subdirectory(${ABACUS_SOURCE_DIR}/module_base ${BASE_BINARY_DIR})
# add parameter
set(PARAMETER_BINARY_DIR "${PROJECT_SOURCE_DIR}/build/parameter")
add_subdirectory(${ABACUS_SOURCE_DIR}/module_parameter ${PARAMETER_BINARY_DIR})
# add orb
set(ORB_BINARY_DIR "${PROJECT_SOURCE_DIR}/build/orb")
add_subdirectory(${ABACUS_SOURCE_DIR}/module_basis/module_ao ${ORB_BINARY_DIR})
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# set RPATH
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# set package name to pyabacus
set(TARGET_PACK pyabacus)
set(CMAKE_INSTALL_RPATH "${PYTHON_SITE_PACKAGES}/${TARGET_PACK}")

# add subdirectories for submodules
add_subdirectory(${PROJECT_SOURCE_DIR}/src/hsolver)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ModuleBase)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ModuleNAO)

