cmake_minimum_required(VERSION 3.15...3.26)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

set(ABACUS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../../source")
set(BASE_PATH "${ABACUS_SOURCE_DIR}/module_base")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake")
# add fftw3
find_package(FFTW3 REQUIRED)
add_compile_definitions(__FFTW3)
find_package(LAPACK REQUIRED)
include_directories(${FFTW3_INCLUDE_DIRS})
message(STATUS "FFTW3_INCLUDE_DIRS: ${FFTW3_INCLUDE_DIRS}")
# add scalapack
find_package(ScaLAPACK REQUIRED)
list(APPEND math_libs FFTW3::FFTW3 LAPACK::LAPACK ScaLAPACK::ScaLAPACK)
include_directories(${BASE_PATH} 
    ${ABACUS_SOURCE_DIR}
    ${ABACUS_SOURCE_DIR}/module_base/module_container)
# add library
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# add base
set(BASE_BINARY_DIR "${PROJECT_SOURCE_DIR}/build/base")
add_subdirectory(${ABACUS_SOURCE_DIR}/module_base ${BASE_BINARY_DIR})
# add orb
set(ORB_BINARY_DIR "${PROJECT_SOURCE_DIR}/build/orb")
set(ENABLE_LCAO ON)
add_subdirectory(${ABACUS_SOURCE_DIR}/module_basis/module_ao ${ORB_BINARY_DIR})
# add nao shared library
file(GLOB _naos ${ABACUS_SOURCE_DIR}/module_basis/module_nao/*.cpp  
${ABACUS_SOURCE_DIR}/module_base/kernels/math_op.cpp
${ABACUS_SOURCE_DIR}/module_psi/kernels/psi_memory_op.cpp
${ABACUS_SOURCE_DIR}/module_io/output_radial.cpp
)
add_library(naopack SHARED 
    ${_naos}
    )
# link math_libs
target_link_libraries(naopack
    base
    container
    orb
    ${math_libs}
)
# list(APPEND _sources ${_naos} ${_bases})
list(APPEND _sources
    # ${ABACUS_SOURCE_DIR}/module_basis/module_nao/radial_collection.cpp
    ${PROJECT_SOURCE_DIR}/src/py_abacus.cpp
    ${PROJECT_SOURCE_DIR}/src/py_base_math.cpp
    ${PROJECT_SOURCE_DIR}/src/py_m_nao.cpp
    )
pybind11_add_module(_core MODULE ${_sources})
target_link_libraries(_core PRIVATE pybind11::headers naopack)
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
install(TARGETS _core naopack DESTINATION pyabacus)

