# Find the oneMKL components used by ABACUS.
#
# ABACUS uses the LP64 BLAS/LAPACK, FFTW3 and, with MPI, BLACS/ScaLAPACK
# interfaces directly.  This module therefore provides complete link closures:
#
#   abacus::mkl              BLAS, LAPACK and FFTW3 compatibility interfaces
#   abacus::mkl_scalapack    abacus::mkl plus BLACS, ScaLAPACK and MPI
#                            (available only when ENABLE_MPI is ON)
#
# Search roots: MKLROOT, MKL_ROOT, or the MKLROOT environment variable.
#
# Optional cache variables:
#   MKL_LINK       AUTO (default), static, or dynamic
#   MKL_MPI        AUTO (default), openmpi, intelmpi, or mpich
#

include(FindPackageHandleStandardArgs)

# Reuse targets configured by this module.  ABACUS deliberately keeps its
# adapter targets separate from the provider-owned MKL:: namespace.
if(TARGET abacus::mkl)
  set(MKL_LIBRARIES abacus::mkl)
  if(ENABLE_MPI)
    if(NOT TARGET abacus::mkl_scalapack)
      message(FATAL_ERROR
        "The existing ABACUS MKL configuration lacks abacus::mkl_scalapack "
        "for an MPI build.")
    endif()
    set(MKL_LIBRARIES abacus::mkl_scalapack)
  endif()
  set(MKL_FOUND TRUE)
  return()
endif()

# ABACUS declares and calls LP64 Fortran-style symbols directly, so ILP64 is not
# ABI-compatible with its integer arguments.
if(DEFINED MKL_INTERFACE)
  string(TOLOWER "${MKL_INTERFACE}" _mkl_integer_interface)
  if(NOT _mkl_integer_interface STREQUAL "lp64")
    message(FATAL_ERROR "ABACUS supports only MKL_INTERFACE=lp64.")
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(_mkl_interface_name mkl_gf_lp64)
else()
  set(_mkl_interface_name mkl_intel_lp64)
endif()

set(_mkl_root_hints "${MKLROOT}" "${MKL_ROOT}" "$ENV{MKLROOT}")
list(REMOVE_ITEM _mkl_root_hints "")
list(REMOVE_DUPLICATES _mkl_root_hints)

set(MKL_LINK AUTO CACHE STRING "oneMKL link mode: AUTO, static, or dynamic")
set_property(CACHE MKL_LINK PROPERTY STRINGS AUTO static dynamic)
string(TOLOWER "${MKL_LINK}" _mkl_link)
if(NOT _mkl_link MATCHES "^(auto|static|dynamic)$")
  message(FATAL_ERROR "MKL_LINK must be AUTO, static, or dynamic.")
endif()

# Keep MKL threading internal: derive it from ABACUS OpenMP support and the
# known compiler/runtime combinations. Unknown OpenMP runtimes use sequential MKL.
set(_mkl_threading sequential)
if(USE_OPENMP)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(_mkl_threading gnu_thread)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(_mkl_threading intel_thread)
  endif()
endif()

if(ENABLE_MPI)
  set(MKL_MPI AUTO CACHE STRING "oneMKL MPI interface: AUTO, openmpi, intelmpi, or mpich")
  set_property(CACHE MKL_MPI PROPERTY STRINGS AUTO openmpi intelmpi mpich)

  string(TOLOWER "${MKL_MPI}" _mkl_mpi)
  if(_mkl_mpi STREQUAL "auto")
    if(CMAKE_CROSSCOMPILING)
      message(
        FATAL_ERROR "CMake cannot auto-determine MKL-BLACS interface correctly with cross compiling. "
                    "Please manually pass -DMKL_MPI=<intelmpi|mpich|openmpi> flag to CMake.")
    endif()
    if("${MPI_CXX_LIBRARY_VERSION_STRING}" MATCHES "Open MPI")
      set(_mkl_mpi openmpi)
    else()
      set(_mkl_mpi intelmpi)
    endif()
  elseif(_mkl_mpi STREQUAL "mpich")
    # MPICH uses Intel-MPI-compatible BLACS on Unix.
    set(_mkl_mpi intelmpi)
  elseif(NOT _mkl_mpi MATCHES "^(openmpi|intelmpi)$")
    message(FATAL_ERROR "MKL_MPI must be AUTO, openmpi, intelmpi, or mpich.")
  endif()

  set(_mkl_blacs_name mkl_blacs_${_mkl_mpi}_lp64)
endif()

# These are result variables of this finder, not user configuration inputs.
foreach(_mkl_result_var IN ITEMS
    MKL_INCLUDE
    MKL_FFTW_INCLUDE
    MKL_INTERFACE_LIB
    MKL_THREAD
    MKL_CORE
    MKL_SCALAPACK
    MKL_BLACS)
  if(DEFINED CACHE{${_mkl_result_var}})
    message(WARNING
      "${_mkl_result_var} is an internal result of FindMKL.cmake and will be ignored.")
    unset(${_mkl_result_var} CACHE)
  endif()
  unset(${_mkl_result_var})
endforeach()

foreach(_mkl_internal_var IN ITEMS
    _abacus_mkl_include
    _abacus_mkl_fftw_include
    _abacus_mkl_interface_lib
    _abacus_mkl_thread
    _abacus_mkl_core
    _abacus_mkl_scalapack
    _abacus_mkl_blacs)
  unset(${_mkl_internal_var})
  unset(${_mkl_internal_var} CACHE)
endforeach()

set(_mkl_saved_suffixes "${CMAKE_FIND_LIBRARY_SUFFIXES}")
if(_mkl_link STREQUAL "static")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
elseif(_mkl_link STREQUAL "dynamic")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
endif()

find_path(_abacus_mkl_include
  NAMES mkl.h
  PATHS ${_mkl_root_hints}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH)

find_path(_abacus_mkl_fftw_include
  NAMES fftw3.h
  PATHS ${_mkl_root_hints}
  PATH_SUFFIXES include/fftw
  NO_DEFAULT_PATH)

find_library(_abacus_mkl_interface_lib
  NAMES ${_mkl_interface_name}
  PATHS ${_mkl_root_hints}
  PATH_SUFFIXES lib/intel64 lib
  NO_DEFAULT_PATH)

find_library(_abacus_mkl_thread
  NAMES mkl_${_mkl_threading}
  PATHS ${_mkl_root_hints}
  PATH_SUFFIXES lib/intel64 lib
  NO_DEFAULT_PATH)

find_library(_abacus_mkl_core
  NAMES mkl_core
  PATHS ${_mkl_root_hints}
  PATH_SUFFIXES lib/intel64 lib
  NO_DEFAULT_PATH)

if(ENABLE_MPI)
  find_library(_abacus_mkl_scalapack
    NAMES mkl_scalapack_lp64
    PATHS ${_mkl_root_hints}
    PATH_SUFFIXES lib/intel64 lib
    NO_DEFAULT_PATH)

  find_library(_abacus_mkl_blacs
    NAMES ${_mkl_blacs_name}
    PATHS ${_mkl_root_hints}
    PATH_SUFFIXES lib/intel64 lib
    NO_DEFAULT_PATH)
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES "${_mkl_saved_suffixes}")

set(MKL_INCLUDE "${_abacus_mkl_include}")
set(MKL_FFTW_INCLUDE "${_abacus_mkl_fftw_include}")
set(MKL_INTERFACE_LIB "${_abacus_mkl_interface_lib}")
set(MKL_THREAD "${_abacus_mkl_thread}")
set(MKL_CORE "${_abacus_mkl_core}")

if(ENABLE_MPI)
  set(MKL_SCALAPACK "${_abacus_mkl_scalapack}")
  set(MKL_BLACS "${_abacus_mkl_blacs}")
endif()

set(_mkl_required_vars MKL_INCLUDE MKL_FFTW_INCLUDE MKL_INTERFACE_LIB MKL_THREAD MKL_CORE)
if(ENABLE_MPI)
  list(APPEND _mkl_required_vars MKL_SCALAPACK MKL_BLACS)
endif()
find_package_handle_standard_args(MKL REQUIRED_VARS ${_mkl_required_vars})

if(MKL_FOUND)
  set(_mkl_libraries ${MKL_INTERFACE_LIB} ${MKL_THREAD} ${MKL_CORE})
  if(ENABLE_MPI)
    message(STATUS "oneMKL BLACS interface: ${_mkl_blacs_name}")
    list(APPEND _mkl_libraries ${MKL_SCALAPACK} ${MKL_BLACS})
  endif()
  set(_mkl_any_static FALSE)
  set(_mkl_all_static TRUE)
  foreach(_mkl_library IN LISTS _mkl_libraries)
    if(_mkl_library MATCHES "\\.a$")
      set(_mkl_any_static TRUE)
    else()
      set(_mkl_all_static FALSE)
    endif()
  endforeach()
  if(_mkl_any_static AND NOT _mkl_all_static)
    message(
      FATAL_ERROR "The selected oneMKL libraries mix static and shared files. "
                  "Choose a consistent MKL_LINK mode or library set.")
  endif()
  if(_mkl_link STREQUAL "static" AND NOT _mkl_all_static)
    message(FATAL_ERROR "MKL_LINK=static did not select static oneMKL libraries.")
  elseif(_mkl_link STREQUAL "dynamic" AND _mkl_all_static)
    message(FATAL_ERROR "MKL_LINK=dynamic did not select shared oneMKL libraries.")
  endif()

  # oneMKL requires pthread even with mkl_sequential
  find_package(Threads REQUIRED)
  set(_mkl_runtime Threads::Threads ${CMAKE_DL_LIBS})
  list(APPEND _mkl_runtime m)
  if(NOT _mkl_threading STREQUAL "sequential")
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    list(APPEND _mkl_runtime OpenMP::OpenMP_CXX)
  endif()

  function(_mkl_link_group output)
    if(_mkl_all_static AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      list(JOIN ARGN "," _mkl_archives)
      set(${output} "-Wl,--start-group,${_mkl_archives},--end-group" PARENT_SCOPE)
    else()
      set(${output} "${ARGN}" PARENT_SCOPE)
    endif()
  endfunction()

  set(_mkl_base_archives ${MKL_INTERFACE_LIB} ${MKL_THREAD} ${MKL_CORE})
  _mkl_link_group(_mkl_base ${_mkl_base_archives})

  add_library(abacus_mkl INTERFACE)
  add_library(abacus::mkl ALIAS abacus_mkl)
  target_include_directories(
    abacus_mkl
    INTERFACE
      "${MKL_INCLUDE}"
      "${MKL_FFTW_INCLUDE}")
  target_link_libraries(abacus_mkl INTERFACE ${_mkl_base} ${_mkl_runtime})

  if(ENABLE_MPI)
    add_library(abacus_mkl_scalapack INTERFACE)
    add_library(abacus::mkl_scalapack ALIAS abacus_mkl_scalapack)
    target_include_directories(
      abacus_mkl_scalapack
      INTERFACE
        "${MKL_INCLUDE}"
        "${MKL_FFTW_INCLUDE}")
    if(_mkl_all_static AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      _mkl_link_group(_mkl_cluster ${MKL_SCALAPACK} ${MKL_BLACS} ${_mkl_base_archives})
      target_link_libraries(
        abacus_mkl_scalapack INTERFACE ${_mkl_cluster} MPI::MPI_CXX ${_mkl_runtime})
    else()
      target_link_libraries(
        abacus_mkl_scalapack INTERFACE
          ${MKL_SCALAPACK} ${MKL_BLACS} abacus::mkl MPI::MPI_CXX)
    endif()
    set(MKL_LIBRARIES abacus::mkl_scalapack)
  else()
    set(MKL_LIBRARIES abacus::mkl)
  endif()
endif()

# Compatibility with packages that consume the standard CMake math targets.
if(TARGET abacus::mkl)
  if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    set_property(TARGET BLAS::BLAS PROPERTY
                 INTERFACE_LINK_LIBRARIES abacus::mkl)
  endif()
  if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_property(TARGET LAPACK::LAPACK PROPERTY
                 INTERFACE_LINK_LIBRARIES abacus::mkl)
  endif()
endif()

mark_as_advanced(
  _abacus_mkl_include
  _abacus_mkl_fftw_include
  _abacus_mkl_interface_lib
  _abacus_mkl_thread
  _abacus_mkl_core
  _abacus_mkl_scalapack
  _abacus_mkl_blacs)
