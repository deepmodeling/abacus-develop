# - Find mkl
# Find the native MKL headers and libraries.
#
#  MKL_INCLUDE_DIRS - where to find mkl.h, etc.
#  MKL_LIBRARIES    - List of libraries when using mkl.
#  MKL_FOUND        - True if mkl found.
#

find_path(MKL_INCLUDE_DIR mkl_service.h HINTS ${MKLROOT}/include)

find_library(MKL_INTEL NAMES mkl_intel_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
find_library(MKL_INTEL_THREAD NAMES mkl_intel_thread HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
find_library(MKL_CORE NAMES mkl_core HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
if(ENABLE_MPI)
  find_library(MKL_SCALAPACK NAMES mkl_scalapack_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
  find_library(MKL_BLACS_INTELMPI NAMES mkl_blacs_intelmpi_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE
# if all listed variables are TRUE

if(ENABLE_MPI)
  find_package_handle_standard_args(IntelMKL DEFAULT_MSG MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_SCALAPACK MKL_BLACS_INTELMPI MKL_INCLUDE_DIR)
else()
  find_package_handle_standard_args(IntelMKL MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_INCLUDE_DIR)
endif()

if(IntelMKL_FOUND)
  if (ENABLE_MPI)
    set(MKL_LIBRARIES ${MKL_INTEL} ${MKL_INTEL_THREAD} ${MKL_CORE} ${MKL_SCALAPACK} ${MKL_BLACS_INTELMPI})
  else()
    set(MKL_LIBRARIES ${MKL_INTEL} ${MKL_INTEL_THREAD} ${MKL_CORE})
  endif()
  set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})

  if(NOT TARGET IntelMKL::INTEL)
    add_library(IntelMKL::INTEL UNKNOWN IMPORTED)
    set_target_properties(IntelMKL::INTEL PROPERTIES
      IMPORTED_LOCATION "${MKL_INTEL}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  if(NOT TARGET IntelMKL::INTEL_THREAD)
    add_library(IntelMKL::INTEL_THREAD UNKNOWN IMPORTED)
    set_target_properties(IntelMKL::INTEL_THREAD PROPERTIES
      IMPORTED_LOCATION "${MKL_INTEL_THREAD}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  if(NOT TARGET IntelMKL::CORE)
    add_library(IntelMKL::CORE UNKNOWN IMPORTED)
    set_target_properties(IntelMKL::CORE PROPERTIES
      IMPORTED_LOCATION "${MKL_CORE}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  if(NOT TARGET IntelMKL::SCALAPACK)
    add_library(IntelMKL::SCALAPACK UNKNOWN IMPORTED)
    set_target_properties(IntelMKL::SCALAPACK PROPERTIES
      IMPORTED_LOCATION "${MKL_SCALAPACK}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  if(NOT TARGET IntelMKL::BLACS_INTELMPI)
    add_library(IntelMKL::BLACS_INTELMPI UNKNOWN IMPORTED)
    set_target_properties(IntelMKL::BLACS_INTELMPI PROPERTIES
      IMPORTED_LOCATION "${MKL_BLACS_INTELMPI}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
  add_library(IntelMKL::MKL INTERFACE IMPORTED)
  if (ENABLE_MPI)
    set_property(TARGET IntelMKL::MKL PROPERTY
    INTERFACE_LINK_LIBRARIES
    "-Wl,--start-group"
    IntelMKL::INTEL IntelMKL::INTEL_THREAD IntelMKL::CORE IntelMKL::SCALAPACK IntelMKL::BLACS_INTELMPI
    "-Wl,--end-group"
    )
  else()
    set_property(TARGET IntelMKL::MKL PROPERTY
    INTERFACE_LINK_LIBRARIES
    "-Wl,--start-group"
    IntelMKL::INTEL IntelMKL::INTEL_THREAD IntelMKL::CORE)
  endif()
endif()

if(ENABLE_MPI)
  mark_as_advanced(MKL_INCLUDE_DIR MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_SCALAPACK MKL_BLACS_INTELMPI)
else()
  mark_as_advanced(MKL_INCLUDE_DIR MKL_INTEL MKL_INTEL_THREAD MKL_CORE)
endif()
