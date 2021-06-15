# - Find FFTW3
# Find the native double precision FFTW3 headers and libraries.
#
#  FFTW3_INCLUDE_DIRS  - Where to find FFTW3 headers.
#  FFTW3_LIBRARIES     - List of libraries when using FFTW3.
#  FFTW3_FOUND         - True if FFTW3 is found.
#

find_path(FFTW3_INCLUDE_DIR fftw3.h)
find_library(FFTW3_LIBRARY NAMES fftw3)

# Handle the QUIET and REQUIRED arguments and
# set FFTW3_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_LIBRARY FFTW3_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(FFTW3_FOUND)
    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})

    if(NOT TARGET FFTW3::FFTW3)
        add_library(FFTW3::FFTW3 UNKNOWN IMPORTED)
        set_target_properties(FFTW3::FFTW3 PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${FFTW3_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
    endif()
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY)

