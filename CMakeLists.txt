########################################
# CMake build system
# This file is part of ABACUS
cmake_minimum_required(VERSION 3.18)
########################################

project(ABACUS
    VERSION 2.2.0
    DESCRIPTION "ABACUS is an electronic structure package based on DFT."
    HOMEPAGE_URL "https://github.com/deepmodeling/abacus-develop"
    LANGUAGES CXX
)

option(ENABLE_DEEPKS "Enable DeePKS functionality" OFF)
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)

if(ENABLE_DEEPKS)
  set(CMAKE_CXX_STANDARD 14)
else()
  set(CMAKE_CXX_STANDARD 11)
endif()
if(ENABLE_ASAN)
  add_compile_options(
    -O1
    -g
    -fsanitize=address
    -fno-omit-frame-pointer
  )
  add_link_options(
    -g
    -fsanitize=address
  )
else()
  add_compile_options(-O2 -g)
endif()

set(ABACUS_BIN_NAME abacus)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)
set(ABACUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(ABACUS_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(ABACUS_BIN_PATH ${CMAKE_CURRENT_BINARY_DIR}/${ABACUS_BIN_NAME})

include_directories(${ABACUS_SOURCE_DIR})

find_package(Boost REQUIRED)
find_package(Cereal REQUIRED)
find_package(ELPA REQUIRED)
find_package(MPI REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

if(DEFINED MKL_DIR)
  set(MKL_DIR_CACHE "${MKL_DIR}")
endif()

include(FetchContent)

include_directories(
  ${Cereal_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${MPI_CXX_INCLUDE_PATH}
)

if(ENABLE_DEEPKS)
  find_package(Torch REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
  FetchContent_Declare(
    libnpy
    GIT_REPOSITORY https://github.com.cnpmjs.org/llohse/libnpy.git
  )
  FetchContent_MakeAvailable(libnpy)
  include_directories(
    ${libnpy_SOURCE_DIR}
    ${TORCH_INCLUDE_DIRS}
  )
  add_definitions(-D__DEEPKS)
endif()

if(DEFINED MKL_DIR_CACHE)
  find_package(IntelMKL REQUIRED)
  add_definitions(-D__MKL -DMKL_ILP64)
  include_directories(${MKL_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS}/fftw)
else()
  find_package(FFTW3 REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(ScaLAPACK REQUIRED)
  include_directories(${FFTW3_INCLUDE_DIRS})
endif()

add_definitions(
    -D__EXX
    -D__FFTW3
    -D__FP
    -D__MPI
    -D__OPENMP
    -D__SELINV
    -D__LCAO
    -DMETIS
    -DEXX_DM=3
    -DEXX_H_COMM=2
    -DTEST_EXX_LCAO=0
    -DTEST_EXX_RADIAL=1
    -DUSE_CEREAL_SERIALIZATION
)

add_subdirectory(source)

add_executable(${ABACUS_BIN_NAME}
    source/main.cpp
)

target_link_libraries(${ABACUS_BIN_NAME}
    base
    cell
    symmetry
    md
    neighbor
    orb
    io
    ions
    lcao
    parallel
    mrrr
    pdiag
    pw
    ri
    driver
    -lm
    ELPA::ELPA
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    Threads::Threads
)

if(DEFINED MKL_DIR_CACHE)
  target_link_libraries(${ABACUS_BIN_NAME}
    -lifcore
    IntelMKL::MKL
    OpenMP::OpenMP_CXX
  )
else()
  target_link_libraries(${ABACUS_BIN_NAME}
    FFTW3::FFTW3
    ScaLAPACK::ScaLAPACK
    LAPACK::LAPACK
    OpenMP::OpenMP_CXX
  )
endif()

# if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
target_link_libraries(${ABACUS_BIN_NAME}
  -lgfortran
)
# endif()

if(ENABLE_DEEPKS)
  target_link_libraries(${ABACUS_BIN_NAME}
    ${TORCH_LIBRARIES}
  )
endif()

install(PROGRAMS ${ABACUS_BIN_PATH}
    TYPE BIN
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)

enable_testing()

add_subdirectory(tests)
