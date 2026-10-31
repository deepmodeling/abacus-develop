cmake_minimum_required(VERSION 3.24)
project(PortableFFTCheck LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(USE_CUFFTDX ON)
set(ABACUS_BIN_NAME check_fused_fft_dispatch)
get_filename_component(ABACUS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../.." ABSOLUTE)
find_package(CUDAToolkit REQUIRED)
find_package(mathdx REQUIRED CONFIG COMPONENTS cufftdx)
if(cufftdx_VERSION VERSION_LESS 1.6)
  message(FATAL_ERROR "The fused FFT checks require cuFFTDx 1.6 or newer.")
endif()

add_compile_definitions(__CUDA=1 __CUFFTDX=1)
include_directories("${ABACUS_SOURCE_DIR}/source"
                    "${ABACUS_SOURCE_DIR}/source/module_base/module_container"
                    ${CUDAToolkit_INCLUDE_DIRS})
add_executable(${ABACUS_BIN_NAME} ../src/check_fused_fft_dispatch.cpp)
target_link_libraries(${ABACUS_BIN_NAME} PRIVATE CUDA::cudart)
add_subdirectory("${ABACUS_SOURCE_DIR}/source/module_basis/module_pw" planewave)
