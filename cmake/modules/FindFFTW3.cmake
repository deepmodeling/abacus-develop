# - Find FFTW3
# Find the native double precision FFTW3 headers and libraries.
#
#  FFTW3_INCLUDE_DIRS  - Where to find FFTW3 headers.
#  FFTW3_LIBRARIES     - List of libraries when using FFTW3.
#  FFTW3_FOUND         - True if FFTW3 is found.
#  FFTW3_VERSION       - Version from the selected library's pkgconfig/fftw3.pc, if available.
#

find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    HINTS ${FFTW3_DIR}
    PATH_SUFFIXES "include"
    )
find_library(FFTW3_LIBRARY
    NAMES fftw3
    HINTS ${FFTW3_DIR}
    PATH_SUFFIXES "lib"
    )

if(ENABLE_FLOAT_FFTW)
  find_library(FFTW3_FLOAT_LIBRARY
      NAMES fftw3f
      HINTS ${FFTW3_DIR}
      PATH_SUFFIXES "lib"
      )
endif()

# Both libfftw3.so and libfftw3_omp.so are required for OpenMP builds.
if (ENABLE_OPENMP)
  find_library(FFTW3_OMP_LIBRARY
      NAMES fftw3_omp
      HINTS ${FFTW3_DIR}
      PATH_SUFFIXES "lib"
      )
endif()

# Handle the QUIET and REQUIRED arguments and
# set FFTW3_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
set(_fftw3_required_vars FFTW3_LIBRARY FFTW3_INCLUDE_DIR)
if(ENABLE_OPENMP)
  list(APPEND _fftw3_required_vars FFTW3_OMP_LIBRARY)
endif()
if(ENABLE_FLOAT_FFTW)
  list(APPEND _fftw3_required_vars FFTW3_FLOAT_LIBRARY)
endif()
find_package_handle_standard_args(FFTW3 DEFAULT_MSG ${_fftw3_required_vars})

# Copy the results to the output variables and target.
set(FFTW3_VERSION "")
if(FFTW3_FOUND)
    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
    if (ENABLE_OPENMP)
        list(APPEND FFTW3_LIBRARIES ${FFTW3_OMP_LIBRARY})
    endif()

    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})

  # Read the literal version from metadata beside the selected library.
  get_filename_component(_fftw3_library_dir "${FFTW3_LIBRARY}" DIRECTORY)
  set(_fftw3_pc "${_fftw3_library_dir}/pkgconfig/fftw3.pc")
  if(EXISTS "${_fftw3_pc}")
    file(STRINGS "${_fftw3_pc}" _fftw3_version_line REGEX "^Version:[ \t]*[0-9]")
    if(_fftw3_version_line)
      string(REGEX REPLACE "^Version:[ \t]*" "" FFTW3_VERSION "${_fftw3_version_line}")
      string(STRIP "${FFTW3_VERSION}" FFTW3_VERSION)
    endif()
  endif()

    if(NOT TARGET FFTW3::FFTW3)
        add_library(FFTW3::FFTW3 UNKNOWN IMPORTED)
        set_target_properties(FFTW3::FFTW3 PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${FFTW3_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
    endif()
    if(ENABLE_FLOAT_FFTW AND NOT TARGET FFTW3::FFTW3_FLOAT)
        add_library(FFTW3::FFTW3_FLOAT UNKNOWN IMPORTED)
        set_target_properties(FFTW3::FFTW3_FLOAT PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                IMPORTED_LOCATION "${FFTW3_FLOAT_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
    endif()
    if (ENABLE_OPENMP)
        if(NOT TARGET FFTW3::FFTW3_OMP)
        add_library(FFTW3::FFTW3_OMP UNKNOWN IMPORTED)
        set_target_properties(FFTW3::FFTW3_OMP PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${FFTW3_OMP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
        endif()
    endif()
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY FFTW3_OMP_LIBRARY)
