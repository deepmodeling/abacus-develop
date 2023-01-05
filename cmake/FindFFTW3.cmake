# - Find FFTW3
# Find the native double precision FFTW3 headers and libraries.
#
#  FFTW3_INCLUDE_DIRS  - Where to find FFTW3 headers.
#  FFTW3_LIBRARIES     - List of libraries when using FFTW3.
#  FFTW3_FOUND         - True if FFTW3 is found.
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
find_library(FFTW3_FLOAT_LIBRARY
    NAMES fftw3f
    HINTS ${FFTW3_DIR}
    PATH_SUFFIXES "lib"
    )

# both libfftw3.so and libfftw3_omp.so should be link in multi-thread term
if (USE_OPENMP)
find_library(FFTW3_OMP_LIBRARY
    NAMES fftw3_omp
    HINTS ${FFTW3_DIR}
    PATH_SUFFIXES "lib"
    )
endif()

# Handle the QUIET and REQUIRED arguments and
# set FFTW3_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
if (USE_OPENMP)
find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_OMP_LIBRARY FFTW3_LIBRARY FFTW3_FLOAT_LIBRARY FFTW3_INCLUDE_DIR)
else()
find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_LIBRARY FFTW3_FLOAT_LIBRARY FFTW3_INCLUDE_DIR)
endif()

# Copy the results to the output variables and target.
if(FFTW3_FOUND)
    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
    if (USE_OPENMP)
        list(APPEND FFTW3_LIBRARIES ${FFTW3_OMP_LIBRARY})
    endif()

    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})

    if(NOT TARGET FFTW3::FFTW3)
        add_library(FFTW3::FFTW3 UNKNOWN IMPORTED)
        set_target_properties(FFTW3::FFTW3 PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${FFTW3_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
    endif()
    if(NOT TARGET FFTW3::FFTW3_FLOAT)
        add_library(FFTW3::FFTW3_FLOAT UNKNOWN IMPORTED)
        set_target_properties(FFTW3::FFTW3_FLOAT PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                IMPORTED_LOCATION "${FFTW3_FLOAT_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
    endif()
    if (USE_OPENMP)
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
