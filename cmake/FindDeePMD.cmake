###############################################################################
# - Find DeePMD
# Find the native DeePMD headers and libraries.
#
#  DeePMD_FOUND        - True if lib is found.
#  DeePMD_LIBRARIES    - List of libraries
#  DeePMD_INCLUDE_DIR  - Where to find DeePMD headers.
#

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
find_library(deepmd_op
    NAMES deepmd_op
    HINTS ${DeePMD_DIR}
    PATH_SUFFIXES "lib"
    )
find_library(deepmd_op_cuda
    NAMES deepmd_op_cuda
    HINTS ${DeePMD_DIR}
    PATH_SUFFIXES "lib"
    )
find_library(tensorflow_cc
    NAMES tensorflow_cc
    HINTS ${DeePMD_DIR}
    PATH_SUFFIXES "lib"
    )

# Handle the QUIET and REQUIRED arguments and
# set DeePMD_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DeePMD DEFAULT_MSG deepmd_cc deepmd_op deepmd_op_cuda tensorflow_cc DeePMD_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(DeePMD_FOUND)
    #set(DeePMD_LIBRARIES ${DeePMD_LIBRARY})
    set(DeePMD_INCLUDE_DIR ${DeePMD_INCLUDE_DIR})

    if(NOT TARGET DeePMD::deepmd_cc)
        add_library(DeePMD::deepmd_cc UNKNOWN IMPORTED)
        set_target_properties(DeePMD::deepmd_cc PROPERTIES
           IMPORTED_LINK_INTERFACE_LANGUAGES "C"
           IMPORTED_LOCATION "${deepmd_cc}"
           INTERFACE_INCLUDE_DIRECTORIES "${DeePMD_INCLUDE_DIR}")
    endif()
    if(NOT TARGET DeePMD::deepmd_op)
        add_library(DeePMD::deepmd_op UNKNOWN IMPORTED)
        set_target_properties(DeePMD::deepmd_op PROPERTIES
           IMPORTED_LINK_INTERFACE_LANGUAGES "C"
           IMPORTED_LOCATION "${deepmd_op}"
           INTERFACE_INCLUDE_DIRECTORIES "${DeePMD_INCLUDE_DIR}")
    endif()
    if(NOT TARGET DeePMD::deepmd_op_cuda)
        add_library(DeePMD::deepmd_op_cuda UNKNOWN IMPORTED)
        set_target_properties(DeePMD::deepmd_op_cuda PROPERTIES
           IMPORTED_LINK_INTERFACE_LANGUAGES "C"
           IMPORTED_LOCATION "${deepmd_op_cuda}"
           INTERFACE_INCLUDE_DIRECTORIES "${DeePMD_INCLUDE_DIR}")
    endif()
    if(NOT TARGET DeePMD::tensorflow_cc)
        add_library(DeePMD::tensorflow_cc UNKNOWN IMPORTED)
        set_target_properties(DeePMD::tensorflow_cc PROPERTIES
           IMPORTED_LINK_INTERFACE_LANGUAGES "C"
           IMPORTED_LOCATION "${tensorflow_cc}"
           INTERFACE_INCLUDE_DIRECTORIES "${DeePMD_INCLUDE_DIR}")
    endif()
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${DeePMD_INCLUDE_DIR})

mark_as_advanced(DeePMD_INCLUDE_DIR deepmd_cc deepmd_op deepmd_op_cuda tensorflow_cc)
