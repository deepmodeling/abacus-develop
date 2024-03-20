###############################################################################
# - Find TensorFlow
# Find the native TensorFlow libraries.
#
#  TensorFlow_FOUND        - True if lib is found.
#  TensorFlow_LIBRARIES    - List of libraries


find_library(tensorflow_cc
    NAMES tensorflow_cc
    HINTS ${TensorFlow_DIR}
    PATH_SUFFIXES "lib"
    )

# Handle the QUIET and REQUIRED arguments and
# set TensorFlow_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TensorFlow DEFAULT_MSG tensorflow_cc)

if(TensorFlow_FOUND)
    if(NOT TARGET TensorFlow::tensorflow_cc)
        add_library(TensorFlow::tensorflow_cc UNKNOWN IMPORTED)
        set_target_properties(TensorFlow::tensorflow_cc PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${tensorflow_cc}")
    endif()
endif()

mark_as_advanced(tensorflow_cc)
