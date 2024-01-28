cmake_minimum_required(VERSION 3.15...3.26)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

set(BASE_PATH "${PROJECT_SOURCE_DIR}/../../source/module_base")
set(ABACUS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../../source")
include_directories(${BASE_PATH} ${ABACUS_SOURCE_DIR})
list(APPEND _sources
    ${ABACUS_SOURCE_DIR}/module_basis/module_nao/numerical_radial.h
    ${ABACUS_SOURCE_DIR}/module_basis/module_nao/numerical_radial.cpp
    ${PROJECT_SOURCE_DIR}/src/py_numerical_radial.cpp)
python_add_library(_core MODULE ${_sources} WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers)
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _core DESTINATION pyabacus)