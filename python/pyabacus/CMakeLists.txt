cmake_minimum_required(VERSION 3.15...3.26)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

set(ABACUS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../../source")
set(BASE_PATH "${ABACUS_SOURCE_DIR}/module_base")
set(NAO_PATH "${ABACUS_SOURCE_DIR}/module_basis/module_nao")
set(ENABLE_LCAO ON)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake")

# add math_libs
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
add_subdirectory(${ABACUS_SOURCE_DIR}/module_basis/module_ao ${ORB_BINARY_DIR})
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# add nao shared library  
list(APPEND _naos
    ${NAO_PATH}/atomic_radials.cpp
    ${NAO_PATH}/beta_radials.cpp
    ${NAO_PATH}/hydrogen_radials.cpp
    ${NAO_PATH}/numerical_radial.cpp
    ${NAO_PATH}/pswfc_radials.cpp
    ${NAO_PATH}/radial_collection.cpp
    ${NAO_PATH}/radial_set.cpp
    ${NAO_PATH}/real_gaunt_table.cpp
    ${NAO_PATH}/sphbes_radials.cpp
    ${NAO_PATH}/two_center_bundle.cpp
    ${NAO_PATH}/two_center_integrator.cpp
    ${NAO_PATH}/two_center_table.cpp
    # dependency
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
    ${PROJECT_SOURCE_DIR}/src/py_abacus.cpp
    ${PROJECT_SOURCE_DIR}/src/py_base_math.cpp
    ${PROJECT_SOURCE_DIR}/src/py_m_nao.cpp
    )
pybind11_add_module(_core MODULE ${_sources})
target_link_libraries(_core PRIVATE pybind11::headers naopack)
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
# set RPATH
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# set package name to pyabacus
set(TARGET_PACK pyabacus)
set(CMAKE_INSTALL_RPATH "${PYTHON_SITE_PACKAGES}/${TARGET_PACK}")
set_target_properties(_core PROPERTIES INSTALL_RPATH "$ORIGIN")
set_target_properties(naopack PROPERTIES INSTALL_RPATH "$ORIGIN")
install(TARGETS _core naopack DESTINATION ${TARGET_PACK})

