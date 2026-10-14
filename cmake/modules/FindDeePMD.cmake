###############################################################################
# - Find DeePMD
# Find the native DeePMD headers and libraries.
#
#  DeePMD_FOUND        - True if lib is found.
#  DeePMDC_FOUND       - True if C API is found.
#  DeePMD_LIBRARIES    - List of libraries
#  DeePMD_INCLUDE_DIR  - Where to find DeePMD headers.
#

# C API
find_path(DeePMD_INCLUDE_C_DIR
    deepmd/deepmd.hpp
    deepmd/c_api.h
    HINTS ${DeePMD_DIR}
    PATH_SUFFIXES "include"
    )
find_library(deepmd_c
    NAMES deepmd_c
    HINTS ${DeePMD_DIR}
    PATH_SUFFIXES "lib"
    )
# C++ API
find_path(DeePMD_INCLUDE_DIR
    deepmd/DeepPot.h
    HINTS ${DeePMD_DIR}
    PATH_SUFFIXES "include"
    )
find_library(deepmd_cc
    NAMES deepmd_cc
    HINTS ${DeePMD_DIR}
    PATH_SUFFIXES "lib"
    )

# Handle the QUIET and REQUIRED arguments and
# set DeePMD_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
if (DeePMD_INCLUDE_C_DIR)
    find_package_handle_standard_args(DeePMD DEFAULT_MSG deepmd_c DeePMD_INCLUDE_C_DIR)
    if (DeePMD_FOUND)
        set(DeePMDC_FOUND TRUE)
        set(DeePMD_INCLUDE_DIR ${DeePMD_INCLUDE_C_DIR})
        if(NOT TARGET DeePMD::deepmd_c)
            add_library(DeePMD::deepmd_c UNKNOWN IMPORTED)
            set_target_properties(DeePMD::deepmd_c PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                IMPORTED_LOCATION "${deepmd_c}"
                INTERFACE_INCLUDE_DIRECTORIES "${DeePMD_INCLUDE_DIR}")
        endif()
    endif()
else()
    find_package_handle_standard_args(DeePMD DEFAULT_MSG deepmd_cc DeePMD_INCLUDE_DIR)
    if(DeePMD_FOUND)
        set(DeePMD_INCLUDE_DIR ${DeePMD_INCLUDE_DIR})
        if(NOT TARGET DeePMD::deepmd_cc)
            add_library(DeePMD::deepmd_cc UNKNOWN IMPORTED)
            set_target_properties(DeePMD::deepmd_cc PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                IMPORTED_LOCATION "${deepmd_cc}"
                INTERFACE_INCLUDE_DIRECTORIES "${DeePMD_INCLUDE_DIR}")
        endif()
    endif()
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${DeePMD_INCLUDE_DIR})

include(CheckCXXSourceCompiles)
if(DeePMD_FOUND)
    # The communication constructor is inline and calls into the DeePMD
    # library, so the feature probe must link against the selected API library.
    set(_abacus_saved_required_libraries ${CMAKE_REQUIRED_LIBRARIES})
    if(DeePMDC_FOUND)
        set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${deepmd_c})
    else()
        set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${deepmd_cc})
    endif()
    if(DeePMDC_FOUND)
        check_cxx_source_compiles(
            "#include <deepmd/deepmd.hpp>
             int main()
             {
                 deepmd::hpp::InputNlist nlist(0, nullptr, nullptr, nullptr,
                                                0, nullptr, nullptr, nullptr,
                                                nullptr, nullptr, nullptr, nullptr, 1);
                 return 0;
             }"
            DeePMD_HAS_COMM_NLIST)
    else()
        check_cxx_source_compiles(
            "#include <deepmd/DeepPot.h>
             int main()
             {
                 deepmd::InputNlist nlist(0, nullptr, nullptr, nullptr,
                                           0, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr, nullptr, 1);
                 return 0;
             }"
            DeePMD_HAS_COMM_NLIST)
    endif()
    set(CMAKE_REQUIRED_LIBRARIES ${_abacus_saved_required_libraries})
    unset(_abacus_saved_required_libraries)
endif()

if(DeePMD_FOUND AND NOT DeePMD_HAS_COMM_NLIST)
    message(WARNING
        "The installed DeePMD API has no communication-capable InputNlist constructor; "
        "message-passing MDCell support requires a newer DeePMD-kit API.")
endif()

mark_as_advanced(DeePMD_INCLUDE_DIR deepmd_c deepmd_cc)
