# - Find mkl
# Find the native MKL headers and libraries.
#
#  MKL_INCLUDE - where to find mkl.h, etc.
#  MKL_FOUND        - True if mkl found.

# find_package(MKL NO_MODULE) # try using official module first
if(NOT TARGET MKL::MKL)

find_path(MKL_INCLUDE mkl_service.h HINTS ${MKLROOT}/include)

find_library(MKL_CORE NAMES mkl_core HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  find_library(MKL_INTERFACE_LIB NAMES mkl_intel_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
  find_library(MKL_THREAD NAMES mkl_intel_thread HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
  find_library(MKL_IOMP5 NAMES iomp5
    HINTS ENV CMPLR_ROOT
    PATH_SUFFIXES lib lib/intel64 linux/compiler/lib/intel64_lin
  )
else()
  find_library(MKL_INTERFACE_LIB NAMES mkl_gf_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
  find_library(MKL_THREAD NAMES mkl_gnu_thread HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
  # With GCC we use system-installed GNU OpenMP
endif()

if(ENABLE_MPI)
  execute_process(COMMAND ${MPI_CXX_COMPILER} --showme:version
                  OUTPUT_VARIABLE MPI_VER_OUT
                  ERROR_VARIABLE MPI_VER_ERR)
  if(MPI_VER_OUT MATCHES "Open MPI" OR MPI_VER_ERR MATCHES "Open MPI")
    set(MKL_BLACS_LIB_NAME "mkl_blacs_openmpi_lp64")
  else()
    set(MKL_BLACS_LIB_NAME "mkl_blacs_intelmpi_lp64")
  endif()
  find_library(MKL_SCALAPACK NAMES mkl_scalapack_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
  find_library(MKL_BLACS NAMES ${MKL_BLACS_LIB_NAME} HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE
# if all listed variables are TRUE

if(ENABLE_MPI)
  find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INTERFACE_LIB MKL_THREAD MKL_CORE MKL_SCALAPACK MKL_BLACS MKL_INCLUDE)
else()
  find_package_handle_standard_args(MKL MKL_INTERFACE_LIB MKL_THREAD MKL_CORE MKL_INCLUDE)
endif()

if(MKL_FOUND)
  if(NOT TARGET MKL::INTERFACE_LIB)
    add_library(MKL::INTERFACE_LIB UNKNOWN IMPORTED)
    set_target_properties(MKL::INTERFACE_LIB PROPERTIES
      IMPORTED_LOCATION "${MKL_INTERFACE_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE}")
  endif()
  if(NOT TARGET MKL::THREAD)
    add_library(MKL::THREAD UNKNOWN IMPORTED)
    set_target_properties(MKL::THREAD PROPERTIES
      IMPORTED_LOCATION "${MKL_THREAD}"
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
  if(ENABLE_MPI AND NOT TARGET MKL::BLACS)
    add_library(MKL::BLACS UNKNOWN IMPORTED)
    set_target_properties(MKL::BLACS PROPERTIES
      IMPORTED_LOCATION "${MKL_BLACS}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE}")
  endif()
  if(MKL_IOMP5 AND NOT TARGET MKL::IOMP5)
    add_library(MKL::IOMP5 UNKNOWN IMPORTED)
    set_target_properties(MKL::IOMP5 PROPERTIES
      IMPORTED_LOCATION "${MKL_IOMP5}")
  endif()
  add_library(MKL::MKL INTERFACE IMPORTED)
  if (ENABLE_MPI)
    set_property(TARGET MKL::MKL PROPERTY
    INTERFACE_LINK_LIBRARIES
    "-Wl,--start-group"
    MKL::INTERFACE_LIB MKL::THREAD MKL::CORE MKL::MKL_SCALAPACK MKL::BLACS
    "-Wl,--end-group"
    )
  else()
    set_property(TARGET MKL::MKL PROPERTY
    INTERFACE_LINK_LIBRARIES
    "-Wl,--start-group"
    MKL::INTERFACE_LIB MKL::THREAD MKL::CORE
    "-Wl,--end-group"
    )
  endif()
  if(TARGET MKL::IOMP5)
    set_property(TARGET MKL::MKL APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES MKL::IOMP5)
  endif()
endif()

if(ENABLE_MPI)
  mark_as_advanced(MKL_INCLUDE MKL_INTERFACE_LIB MKL_THREAD MKL_CORE MKL_SCALAPACK MKL_BLACS)
else()
  mark_as_advanced(MKL_INCLUDE MKL_INTERFACE_LIB MKL_THREAD MKL_CORE)
endif()

endif() # MKL::MKL



# In oneAPI 2022, MKL_SCALAPACK might not be linked properly
if(NOT TARGET MKL::MKL_SCALAPACK)
  find_library(MKL_SCALAPACK NAMES mkl_scalapack_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
  message(STATUS "Found MKL_SCALAPACK: ${MKL_SCALAPACK}")
  if(MKL_SCALAPACK)
    # create an IMPORTED target that points to the discovered library file
    add_library(MKL::MKL_SCALAPACK UNKNOWN IMPORTED)
    set_target_properties(MKL::MKL_SCALAPACK PROPERTIES
      IMPORTED_LOCATION "${MKL_SCALAPACK}"
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE}"
    )
  endif()
endif()
