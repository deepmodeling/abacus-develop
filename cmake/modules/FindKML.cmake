# - Find Huawei Kunpeng Math Library (KML)
#
# This module finds the KML linear-algebra libraries. KML installs optimized
# variants below architecture- and threading-specific directories, so those
# choices are made here rather than at each consumer.
#
# Imported targets:
#   KML::BLAS         KML BLAS
#   KML::LAPACK       KML LAPACK, with KML::BLAS transitively linked
#   KML::ScaLAPACK    KML ScaLAPACK, with KML::LAPACK transitively linked
#   KML::FFTW3        KML FFTW-compatible double-precision FFT interface
#   KML::FFTW3_FLOAT  KML FFTW-compatible single-precision FFT interface
#
# If libkml_rt is present, KML::BLAS also propagates it as a runtime
# dependency.
#
# Cache variables:
#   KML_ROOT            KML installation prefix
#   KML_ARCH            KML library variant: neon, sve, or sve512
#   KML_BLAS_THREADING  kblas variant: auto, multi, locking, or nolocking
#                        (default: auto)
#
# The default threading selection uses the caller's USE_OPENMP option when it
# is available: multi for OpenMP builds and nolocking otherwise. Projects can
# select a KML_BLAS_THREADING variant explicitly.
#
# KML_ARCH=sve512 uses the lib/sme KBLAS directory; the remaining KML
# libraries use lib/sve512.

include(FindPackageHandleStandardArgs)

set(KML_ROOT "" CACHE PATH "KML installation prefix")

set(KML_ARCH "neon" CACHE STRING "KML library variant (neon, sve, or sve512)")
set_property(CACHE KML_ARCH PROPERTY STRINGS neon sve sve512)

set(KML_BLAS_THREADING "auto" CACHE STRING
    "KML kblas variant (auto, multi, locking, or nolocking)")
set_property(CACHE KML_BLAS_THREADING PROPERTY STRINGS
             auto multi locking nolocking)

set(_kml_arch_variants neon sve sve512)
if(NOT KML_ARCH IN_LIST _kml_arch_variants)
  message(FATAL_ERROR "KML_ARCH must be one of: ${_kml_arch_variants}")
endif()

set(_kml_thread_variants multi locking nolocking)
if(KML_BLAS_THREADING STREQUAL "auto")
  if(DEFINED USE_OPENMP AND USE_OPENMP)
    set(_kml_blas_threading multi)
  else()
    set(_kml_blas_threading nolocking)
  endif()
else()
  set(_kml_blas_threading "${KML_BLAS_THREADING}")
endif()

if(NOT _kml_blas_threading IN_LIST _kml_thread_variants)
  message(FATAL_ERROR
    "KML_BLAS_THREADING must be auto or one of: ${_kml_thread_variants}")
endif()

set(_kml_blas_arch "${KML_ARCH}")
if(KML_ARCH STREQUAL "sve512")
  set(_kml_blas_arch sme)
endif()

set(_kml_prefix_hints)
if(KML_ROOT)
  list(APPEND _kml_prefix_hints "${KML_ROOT}")
endif()
if(DEFINED ENV{KML_ROOT})
  list(APPEND _kml_prefix_hints "$ENV{KML_ROOT}")
endif()
list(APPEND _kml_prefix_hints /usr/local/kml)

# Check if an explicitly selected compiler-specific KML prefix matches the compiler.
set(_kml_explicit_root "${KML_ROOT}")
if(NOT _kml_explicit_root AND DEFINED ENV{KML_ROOT})
  set(_kml_explicit_root "$ENV{KML_ROOT}")
endif()
get_filename_component(_kml_root_name "${_kml_explicit_root}" NAME)
if(_kml_root_name STREQUAL "gcc")
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
      "KML_ROOT points to the GCC KML bundle, but the C++ compiler is "
      "${CMAKE_CXX_COMPILER_ID}.")
  endif()
elseif(_kml_root_name STREQUAL "bisheng")
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR
      "KML_ROOT points to the BiShengLLVM KML bundle, but the C++ compiler is "
      "${CMAKE_CXX_COMPILER_ID}.")
  endif()
endif()

find_path(KML_INCLUDE_DIR
  NAMES kblas.h klapack.h kscalapack.h
  HINTS ${_kml_prefix_hints}
  PATH_SUFFIXES include)

if(KML_INCLUDE_DIR)
  get_filename_component(_kml_prefix "${KML_INCLUDE_DIR}" DIRECTORY)

  find_library(KML_RUNTIME_LIBRARY
    NAMES kml_rt
    HINTS "${_kml_prefix}"
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH)

  find_library(KML_BLAS_LIBRARY
    NAMES kblas
    HINTS "${_kml_prefix}"
    PATH_SUFFIXES "lib/${_kml_blas_arch}/kblas/${_kml_blas_threading}"
    NO_DEFAULT_PATH)

  find_library(KML_LAPACK_LIBRARY
    NAMES klapack_full
    HINTS "${_kml_prefix}"
    PATH_SUFFIXES "lib/${KML_ARCH}"
    NO_DEFAULT_PATH)

  find_library(KML_SCALAPACK_LIBRARY
    NAMES kscalapack_full
    HINTS "${_kml_prefix}"
    PATH_SUFFIXES "lib/${KML_ARCH}"
    NO_DEFAULT_PATH)

  find_library(KML_FFTW3_LIBRARY
    NAMES fftw3
    HINTS "${_kml_prefix}"
    PATH_SUFFIXES lib/noarch
    NO_DEFAULT_PATH)

  find_library(KML_KFFT_LIBRARY
    NAMES kfft
    HINTS "${_kml_prefix}"
    PATH_SUFFIXES "lib/${KML_ARCH}"
    NO_DEFAULT_PATH)

  find_library(KML_FFTW3_FLOAT_LIBRARY
    NAMES fftw3f
    HINTS "${_kml_prefix}"
    PATH_SUFFIXES lib/noarch
    NO_DEFAULT_PATH)

  find_library(KML_KFFTF_LIBRARY
    NAMES kfftf
    HINTS "${_kml_prefix}"
    PATH_SUFFIXES "lib/${KML_ARCH}"
    NO_DEFAULT_PATH)
endif()

set(KML_BLAS_FOUND FALSE)
if(KML_INCLUDE_DIR AND KML_BLAS_LIBRARY)
  set(KML_BLAS_FOUND TRUE)
endif()

set(KML_LAPACK_FOUND FALSE)
if(KML_BLAS_FOUND AND KML_LAPACK_LIBRARY)
  set(KML_LAPACK_FOUND TRUE)
endif()

set(KML_ScaLAPACK_FOUND FALSE)
if(KML_LAPACK_FOUND AND KML_SCALAPACK_LIBRARY)
  set(KML_ScaLAPACK_FOUND TRUE)
endif()
set(KML_SCALAPACK_FOUND "${KML_ScaLAPACK_FOUND}")

set(KML_FFTW3_FOUND FALSE)
if(KML_INCLUDE_DIR AND EXISTS "${KML_INCLUDE_DIR}/fftw3.h" AND
   KML_FFTW3_LIBRARY AND KML_KFFT_LIBRARY)
  set(KML_FFTW3_FOUND TRUE)
endif()

set(KML_FFTW3_FLOAT_FOUND FALSE)
if(KML_INCLUDE_DIR AND EXISTS "${KML_INCLUDE_DIR}/fftw3.h" AND
   KML_FFTW3_FLOAT_LIBRARY AND KML_KFFTF_LIBRARY)
  set(KML_FFTW3_FLOAT_FOUND TRUE)
endif()

set(_kml_required_vars KML_INCLUDE_DIR)
if(KML_FIND_COMPONENTS)
  foreach(_kml_component IN LISTS KML_FIND_COMPONENTS)
    if(_kml_component STREQUAL "BLAS")
      list(APPEND _kml_required_vars KML_BLAS_LIBRARY)
    elseif(_kml_component STREQUAL "LAPACK")
      list(APPEND _kml_required_vars KML_LAPACK_LIBRARY KML_BLAS_LIBRARY)
    elseif(_kml_component STREQUAL "ScaLAPACK" OR _kml_component STREQUAL "SCALAPACK")
      set(KML_${_kml_component}_FOUND "${KML_ScaLAPACK_FOUND}")
      list(APPEND _kml_required_vars
           KML_SCALAPACK_LIBRARY KML_LAPACK_LIBRARY KML_BLAS_LIBRARY)
    elseif(_kml_component STREQUAL "FFTW3")
      list(APPEND _kml_required_vars KML_FFTW3_LIBRARY KML_KFFT_LIBRARY)
    elseif(_kml_component STREQUAL "FFTW3_FLOAT")
      list(APPEND _kml_required_vars
           KML_FFTW3_FLOAT_LIBRARY KML_KFFTF_LIBRARY)
    else()
      set(KML_${_kml_component}_FOUND FALSE)
    endif()
  endforeach()
else()
  list(APPEND _kml_required_vars KML_BLAS_LIBRARY KML_LAPACK_LIBRARY)
endif()
list(REMOVE_DUPLICATES _kml_required_vars)

find_package_handle_standard_args(KML
  REQUIRED_VARS ${_kml_required_vars}
  HANDLE_COMPONENTS)

if(KML_FOUND)
  set(KML_INCLUDE_DIRS "${KML_INCLUDE_DIR}")

  if(KML_RUNTIME_LIBRARY AND NOT TARGET KML::Runtime)
    add_library(KML::Runtime UNKNOWN IMPORTED)
    set_target_properties(KML::Runtime PROPERTIES
      IMPORTED_LOCATION "${KML_RUNTIME_LIBRARY}")
  endif()

  if(KML_BLAS_FOUND AND NOT TARGET KML::BLAS)
    add_library(KML::BLAS UNKNOWN IMPORTED)
    set_target_properties(KML::BLAS PROPERTIES
      IMPORTED_LOCATION "${KML_BLAS_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${KML_INCLUDE_DIR}")
    if(TARGET KML::Runtime)
      set_property(TARGET KML::BLAS APPEND PROPERTY
                   INTERFACE_LINK_LIBRARIES KML::Runtime)
    endif()
  endif()

  if(KML_LAPACK_FOUND AND NOT TARGET KML::LAPACK)
    add_library(KML::LAPACK UNKNOWN IMPORTED)
    set_target_properties(KML::LAPACK PROPERTIES
      IMPORTED_LOCATION "${KML_LAPACK_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${KML_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES KML::BLAS)
  endif()

  if(KML_ScaLAPACK_FOUND AND NOT TARGET KML::ScaLAPACK)
    add_library(KML::ScaLAPACK UNKNOWN IMPORTED)
    set_target_properties(KML::ScaLAPACK PROPERTIES
      IMPORTED_LOCATION "${KML_SCALAPACK_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${KML_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES KML::LAPACK)
  endif()

  if(KML_FFTW3_FOUND AND NOT TARGET KML::FFTW3)
    add_library(KML::FFTW3 UNKNOWN IMPORTED)
    set_target_properties(KML::FFTW3 PROPERTIES
      IMPORTED_LOCATION "${KML_FFTW3_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${KML_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${KML_KFFT_LIBRARY}")
  endif()

  if(KML_FFTW3_FLOAT_FOUND AND NOT TARGET KML::FFTW3_FLOAT)
    add_library(KML::FFTW3_FLOAT UNKNOWN IMPORTED)
    set_target_properties(KML::FFTW3_FLOAT PROPERTIES
      IMPORTED_LOCATION "${KML_FFTW3_FLOAT_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${KML_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${KML_KFFTF_LIBRARY}")
  endif()

  if(TARGET KML::LAPACK)
    set(KML_LIBRARIES KML::LAPACK)
  elseif(TARGET KML::BLAS)
    set(KML_LIBRARIES KML::BLAS)
  endif()
endif()

# Compatibility with packages that consume the standard CMake math targets.
if(TARGET KML::BLAS AND NOT TARGET BLAS::BLAS)
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  set_property(TARGET BLAS::BLAS PROPERTY
               INTERFACE_LINK_LIBRARIES KML::BLAS)
endif()
if(TARGET KML::LAPACK AND NOT TARGET LAPACK::LAPACK)
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_property(TARGET LAPACK::LAPACK PROPERTY
               INTERFACE_LINK_LIBRARIES KML::LAPACK)
endif()

mark_as_advanced(
  KML_INCLUDE_DIR
  KML_RUNTIME_LIBRARY
  KML_BLAS_LIBRARY
  KML_LAPACK_LIBRARY
  KML_SCALAPACK_LIBRARY
  KML_FFTW3_LIBRARY
  KML_KFFT_LIBRARY
  KML_FFTW3_FLOAT_LIBRARY
  KML_KFFTF_LIBRARY)
