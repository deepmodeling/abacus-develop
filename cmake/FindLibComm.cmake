###############################################################################
# - Find LibComm
# Find the native LibComm files.
#
#  LIBCOMM_FOUND - True if LibComm is found.
#  LIBCOMM_DIR - Where to find LibComm files.

find_path(LIBCOMM_DIR
    include/Comm/Comm_Tools.h
    HINTS ${LIBCOMM_DIR}
    HINTS ${LibComm_DIR}
    HINTS ${libcomm_DIR}
)

if(NOT LIBCOMM_DIR)
    include(FetchContent)
    FetchContent_Declare(
        LibComm
        URL https://codeload.github.com/abacusmodeling/LibComm/tar.gz/965bf90713ce5598bae4fd957d417065f77ceb42
    )
    FetchContent_Populate(LibComm)
    set(LIBCOMM_DIR ${libcomm_SOURCE_DIR})
endif()
# Handle the QUIET and REQUIRED arguments and
# set LIBCOMM_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibComm DEFAULT_MSG LIBCOMM_DIR)

# Copy the results to the output variables and target.
mark_as_advanced(LIBCOMM_DIR)
