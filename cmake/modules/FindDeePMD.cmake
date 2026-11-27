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
find_path(
  DeePMD_INCLUDE_C_DIR deepmd/deepmd.hpp deepmd/c_api.h
  HINTS ${DeePMD_DIR}
  PATH_SUFFIXES "include")
find_library(
  deepmd_c
  NAMES deepmd_c
  HINTS ${DeePMD_DIR}
  PATH_SUFFIXES "lib")
# C++ API
find_path(
  DeePMD_INCLUDE_DIR deepmd/DeepPot.h
  HINTS ${DeePMD_DIR}
  PATH_SUFFIXES "include")
find_library(
  deepmd_cc
  NAMES deepmd_cc
  HINTS ${DeePMD_DIR}
  PATH_SUFFIXES "lib")

# Handle the QUIET and REQUIRED arguments and
# set DeePMD_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
include(CheckCXXSourceCompiles)
if(DeePMD_INCLUDE_DIR AND deepmd_cc)
  set(_abacus_saved_required_includes ${CMAKE_REQUIRED_INCLUDES})
  set(_abacus_saved_required_libraries ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${DeePMD_INCLUDE_DIR})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${deepmd_cc})
  check_cxx_source_compiles(
    "#include <deepmd/DeepPot.h>
     int main() { deepmd::DeepPot model; return model.is_aparam_nall(); }"
    DeePMD_CXX_API_USABLE)
  set(CMAKE_REQUIRED_INCLUDES ${_abacus_saved_required_includes})
  set(CMAKE_REQUIRED_LIBRARIES ${_abacus_saved_required_libraries})
  unset(_abacus_saved_required_includes)
  unset(_abacus_saved_required_libraries)
endif()
# The C++ API exposes is_aparam_nall() on an initialized model. Prefer it
# when both libraries are installed to avoid loading a second model for this
# property in the C API compatibility path.
if(DeePMD_CXX_API_USABLE)
  find_package_handle_standard_args(DeePMD DEFAULT_MSG deepmd_cc
                                    DeePMD_INCLUDE_DIR)
  if(DeePMD_FOUND)
    set(DeePMDC_FOUND FALSE)
    if(NOT TARGET DeePMD::deepmd_cc)
      add_library(DeePMD::deepmd_cc UNKNOWN IMPORTED)
      set_target_properties(
        DeePMD::deepmd_cc
        PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                   IMPORTED_LOCATION "${deepmd_cc}"
                   INTERFACE_INCLUDE_DIRECTORIES "${DeePMD_INCLUDE_DIR}")
    endif()
  endif()
else()
  find_package_handle_standard_args(DeePMD DEFAULT_MSG deepmd_c
                                    DeePMD_INCLUDE_C_DIR)
  if(DeePMD_FOUND)
    set(DeePMDC_FOUND TRUE)
    set(DeePMD_INCLUDE_DIR ${DeePMD_INCLUDE_C_DIR})
    if(NOT TARGET DeePMD::deepmd_c)
      add_library(DeePMD::deepmd_c UNKNOWN IMPORTED)
      set_target_properties(
        DeePMD::deepmd_c
        PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                   IMPORTED_LOCATION "${deepmd_c}"
                   INTERFACE_INCLUDE_DIRECTORIES "${DeePMD_INCLUDE_DIR}")
    endif()
  endif()
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${DeePMD_INCLUDE_DIR})

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
    set(_abacus_dp_nlist_type "deepmd::hpp::InputNlist")
    set(_abacus_dp_nlist_header "#include <deepmd/deepmd.hpp>")
  else()
    set(_abacus_dp_nlist_type "deepmd::InputNlist")
    set(_abacus_dp_nlist_header "#include <deepmd/DeepPot.h>")
  endif()
  check_cxx_source_compiles(
    "${_abacus_dp_nlist_header}
         int main()
         {
             ${_abacus_dp_nlist_type} nlist(0, nullptr, nullptr, nullptr,
                                            0, nullptr, nullptr, nullptr,
                                            nullptr, nullptr, nullptr, nullptr, 1);
             return 0;
         }"
    DeePMD_HAS_COMM_NLIST_13)
  if(DeePMD_HAS_COMM_NLIST_13)
    set(DeePMD_INPUT_NLIST_COMM_ARGS 13)
  else()
    check_cxx_source_compiles(
      "${_abacus_dp_nlist_header}
             int main()
             {
                 ${_abacus_dp_nlist_type} nlist(0, nullptr, nullptr, nullptr,
                                                0, nullptr, nullptr, nullptr,
                                                nullptr, nullptr, nullptr, nullptr);
                 return 0;
             }"
      DeePMD_HAS_COMM_NLIST_12)
    if(DeePMD_HAS_COMM_NLIST_12)
      set(DeePMD_INPUT_NLIST_COMM_ARGS 12)
    else()
      set(DeePMD_INPUT_NLIST_COMM_ARGS 0)
    endif()
  endif()
  check_cxx_source_compiles(
    "${_abacus_dp_nlist_header}
         int main()
         {
             ${_abacus_dp_nlist_type} nlist(0, nullptr, nullptr, nullptr);
             nlist.set_mapping(nullptr);
             return 0;
         }"
    DeePMD_HAS_NLIST_MAPPING)

  # The InputNlist constructor can exist even when the PyTorch operator was
  # built without USE_MPI. In that build its embedding exchange is local only.
  # Require a direct MPI dependency of the operator before enabling the
  # distributed inference path.
  unset(DeePMD_PT_OPERATOR CACHE)
  find_file(
    DeePMD_PT_OPERATOR
    NAMES libdeepmd_op_pt.so
    HINTS ${DeePMD_DIR}
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH)
  set(DeePMD_PT_OPERATOR_MPI FALSE)
  if(DeePMD_PT_OPERATOR AND CMAKE_READELF)
    execute_process(
      COMMAND ${CMAKE_READELF} -d ${DeePMD_PT_OPERATOR}
      RESULT_VARIABLE _abacus_dp_readelf_result
      OUTPUT_VARIABLE _abacus_dp_operator_dynamic)
    if(_abacus_dp_readelf_result EQUAL 0
       AND _abacus_dp_operator_dynamic MATCHES "NEEDED[^\n]*libmpi\\.so")
      set(DeePMD_PT_OPERATOR_MPI TRUE)
    endif()
    unset(_abacus_dp_readelf_result)
    unset(_abacus_dp_operator_dynamic)
  endif()
  if(NOT DeePMD_PT_OPERATOR_MPI)
    set(DeePMD_INPUT_NLIST_COMM_ARGS 0)
  endif()
  set(CMAKE_REQUIRED_LIBRARIES ${_abacus_saved_required_libraries})
  unset(_abacus_saved_required_libraries)
  unset(_abacus_dp_nlist_type)
  unset(_abacus_dp_nlist_header)
endif()

if(DeePMD_FOUND AND NOT DeePMD_PT_OPERATOR_MPI)
  message(WARNING
          "Distributed PyTorch DeePMD inference is disabled: the operator's "
          "direct libmpi.so dependency could not be verified.")
elseif(DeePMD_FOUND
       AND DEFINED DeePMD_INPUT_NLIST_COMM_ARGS
       AND DeePMD_INPUT_NLIST_COMM_ARGS EQUAL 0)
  message(WARNING
          "Distributed PyTorch DeePMD inference is disabled: communication-capable "
          "InputNlist is unavailable.")
endif()

mark_as_advanced(DeePMD_INCLUDE_DIR deepmd_c deepmd_cc)
