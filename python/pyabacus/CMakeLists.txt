cmake_minimum_required(VERSION 3.15...3.26)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
# add fftw3
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake")
find_package(FFTW3 REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(${FFTW3_INCLUDE_DIRS})
list(APPEND math_libs FFTW3::FFTW3 LAPACK::LAPACK)

set(ABACUS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../../source")
include_directories(${BASE_PATH} ${ABACUS_SOURCE_DIR})
list(APPEND _sources
    #${ABACUS_SOURCE_DIR}/module_basis/module_nao/numerical_radial.h
    #${ABACUS_SOURCE_DIR}/module_basis/module_nao/numerical_radial.cpp
    ${ABACUS_SOURCE_DIR}/module_base/constants.h
    ${ABACUS_SOURCE_DIR}/module_base/math_sphbes.h
    ${ABACUS_SOURCE_DIR}/module_base/math_sphbes.cpp
    ${ABACUS_SOURCE_DIR}/module_base/math_integral.h
    ${ABACUS_SOURCE_DIR}/module_base/math_integral.cpp
    ${ABACUS_SOURCE_DIR}/module_base/spherical_bessel_transformer.h
    ${ABACUS_SOURCE_DIR}/module_base/spherical_bessel_transformer.cpp
    ${PROJECT_SOURCE_DIR}/src/py_abacus.cpp
    #${PROJECT_SOURCE_DIR}/src/py_numerical_radial.cpp
    ${PROJECT_SOURCE_DIR}/src/py_math_base.cpp
    )
python_add_library(_core MODULE ${_sources} WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers ${math_libs})
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _core DESTINATION pyabacus)