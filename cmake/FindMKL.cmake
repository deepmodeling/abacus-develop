# - Find mkl
# Find the native MKL headers and libraries.
#
#  MKL_INCLUDE - where to find mkl.h, etc.
#  MKL_FOUND        - True if mkl found.

find_package(MKL NO_MODULE) # try using official module first
if(NOT TARGET MKL::MKL)

find_path(MKL_INCLUDE mkl_service.h HINTS ${MKLROOT}/include)

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
  find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_SCALAPACK MKL_BLACS_INTELMPI MKL_INCLUDE)
else()
  find_package_handle_standard_args(MKL MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_INCLUDE)
endif()

if(MKL_FOUND)
  if(NOT TARGET MKL::INTEL)
    add_library(MKL::INTEL UNKNOWN IMPORTED)
    set_target_properties(MKL::INTEL PROPERTIES
      IMPORTED_LOCATION "${MKL_INTEL}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE}")
  endif()
  if(NOT TARGET MKL::INTEL_THREAD)
    add_library(MKL::INTEL_THREAD UNKNOWN IMPORTED)
    set_target_properties(MKL::INTEL_THREAD PROPERTIES
      IMPORTED_LOCATION "${MKL_INTEL_THREAD}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE}")
  endif()
  if(NOT TARGET MKL::CORE)
    add_library(MKL::CORE UNKNOWN IMPORTED)
    set_target_properties(MKL::CORE PROPERTIES
      IMPORTED_LOCATION "${MKL_CORE}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE}")
  endif()
  if(NOT TARGET MKL::MKL_SCALAPACK)
    add_library(MKL::MKL_SCALAPACK UNKNOWN IMPORTED)
    set_target_properties(MKL::MKL_SCALAPACK PROPERTIES
      IMPORTED_LOCATION "${MKL_SCALAPACK}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE}")
  endif()
  if(NOT TARGET MKL::BLACS_INTELMPI)
    add_library(MKL::BLACS_INTELMPI UNKNOWN IMPORTED)
    set_target_properties(MKL::BLACS_INTELMPI PROPERTIES
      IMPORTED_LOCATION "${MKL_BLACS_INTELMPI}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE}")
  endif()
  add_library(MKL::MKL INTERFACE IMPORTED)
  if (ENABLE_MPI)
    set_property(TARGET MKL::MKL PROPERTY
    INTERFACE_LINK_LIBRARIES
    "-Wl,--start-group"
    MKL::INTEL MKL::INTEL_THREAD MKL::CORE MKL::MKL_SCALAPACK MKL::BLACS_INTELMPI
    "-Wl,--end-group"
    )
  else()
    set_property(TARGET MKL::MKL PROPERTY
    INTERFACE_LINK_LIBRARIES
    "-Wl,--start-group"
    MKL::INTEL MKL::INTEL_THREAD MKL::CORE)
  endif()
endif()

if(ENABLE_MPI)
  mark_as_advanced(MKL_INCLUDE MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_SCALAPACK MKL_BLACS_INTELMPI)
else()
  mark_as_advanced(MKL_INCLUDE MKL_INTEL MKL_INTEL_THREAD MKL_CORE)
endif()

endif() # MKL::MKL

# For compatibility with legacy libpaw_interface CMakeLists.txt
if(TARGET MKL::MKL)
  add_library(IntelMKL::MKL ALIAS MKL::MKL)
endif()

# In oneAPI 2022, MKL_SCALAPACK might not be linked properly
if(NOT TARGET MKL::MKL_SCALAPACK)
  find_library(MKL_SCALAPACK NAMES mkl_scalapack_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
  message(STATUS "Found MKL_SCALAPACK: ${MKL_SCALAPACK}")
  add_library(MKL::MKL_SCALAPACK OBJECT IMPORTED MKL_SCALAPACK)
endif()
