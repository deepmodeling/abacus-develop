cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(nnof)

#list(APPEND CMAKE_PREFIX_PATH /path/to/Torch)
#list(APPEND libnpy_INCLUDE_DIR /path/to/libnpy)

find_package(Torch REQUIRED)
include_directories(${TORCH_INCLUDE_DIRS})

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
    include_directories(${libnpy_INCLUDE_DIR})
endif()
include_directories(${libnpy_SOURCE_DIR}/include)

add_executable(nnof main.cpp data.cpp nn_of.cpp grid.cpp input.cpp kernel.cpp pauli_potential.cpp train_kedf.cpp)
target_link_libraries(nnof "${TORCH_LIBRARIES}")
set_property(TARGET nnof PROPERTY CXX_STANDARD 14)

#cmake -B build -Dlibnpy_INCLUDE_DIR=/path/to/libnpy/ -DTorch_DIR=/path/to/Torch/
