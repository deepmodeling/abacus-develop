# Find the native DeePMD headers and API libraries. DeePMD_DIR is an
# installation prefix; CMAKE_PREFIX_PATH and DeePMD_ROOT are also supported.
# ABACUS_DEEPMD_TARGET names the selected C++ API or C-wrapper target.
find_path(
  DeePMD_INCLUDE_DIR deepmd/DeepPot.h
  HINTS ${DeePMD_DIR} ENV DeePMD_DIR
  PATH_SUFFIXES include)
find_path(
  DeePMD_INCLUDE_C_DIR deepmd/deepmd.hpp
  HINTS ${DeePMD_DIR} ENV DeePMD_DIR
  PATH_SUFFIXES include)
find_library(
  deepmd_cc
  NAMES deepmd_cc
  HINTS ${DeePMD_DIR} ENV DeePMD_DIR
  PATH_SUFFIXES lib lib64)
find_library(
  deepmd_c
  NAMES deepmd_c
  HINTS ${DeePMD_DIR} ENV DeePMD_DIR
  PATH_SUFFIXES lib lib64)

if(deepmd_cc
   AND DeePMD_INCLUDE_DIR
   AND NOT TARGET DeePMD::deepmd_cc)
  add_library(DeePMD::deepmd_cc UNKNOWN IMPORTED)
  set_target_properties(
    DeePMD::deepmd_cc
    PROPERTIES IMPORTED_LOCATION "${deepmd_cc}" INTERFACE_INCLUDE_DIRECTORIES
                                                "${DeePMD_INCLUDE_DIR}")
endif()
if(deepmd_c
   AND DeePMD_INCLUDE_C_DIR
   AND NOT TARGET DeePMD::deepmd_c)
  add_library(DeePMD::deepmd_c UNKNOWN IMPORTED)
  set_target_properties(
    DeePMD::deepmd_c
    PROPERTIES IMPORTED_LOCATION "${deepmd_c}" INTERFACE_INCLUDE_DIRECTORIES
                                               "${DeePMD_INCLUDE_C_DIR}")
endif()

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state(RESET)

# Do not reuse results from a different DeePMD installation or backend.
foreach(_abacus_dp_check IN
        ITEMS DeePMD_CXX_API_USABLE DeePMD_HAS_COMM_NLIST_13
              DeePMD_HAS_COMM_NLIST_12 DeePMD_HAS_NLIST_MAPPING)
  unset(${_abacus_dp_check} CACHE)
endforeach()
if(_enable_torch)
  set(CMAKE_REQUIRED_FLAGS "${TORCH_CXX_FLAGS}")
endif()

if(TARGET DeePMD::deepmd_cc)
  set(CMAKE_REQUIRED_LIBRARIES DeePMD::deepmd_cc)
  check_cxx_source_compiles(
    "#include <deepmd/DeepPot.h>
     int main() { deepmd::DeepPot model; return model.is_aparam_nall(); }"
    DeePMD_CXX_API_USABLE)
endif()

# Prefer the C++ API when it exposes the required property, retaining the
# C-wrapper fallback for installations with an older C++ interface.
set(ABACUS_DEEPMD_TARGET)
set(DeePMDC_FOUND FALSE)
if(DeePMD_CXX_API_USABLE)
  set(ABACUS_DEEPMD_TARGET DeePMD::deepmd_cc)
  set(_abacus_dp_nlist_type "deepmd::InputNlist")
  set(_abacus_dp_nlist_header "#include <deepmd/DeepPot.h>")
elseif(TARGET DeePMD::deepmd_c)
  set(ABACUS_DEEPMD_TARGET DeePMD::deepmd_c)
  set(DeePMDC_FOUND TRUE)
  set(_abacus_dp_nlist_type "deepmd::hpp::InputNlist")
  set(_abacus_dp_nlist_header "#include <deepmd/deepmd.hpp>")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  DeePMD REQUIRED_VARS ABACUS_DEEPMD_TARGET REASON_FAILURE_MESSAGE
                       "A usable C++ API or the C API is required.")

if(DeePMD_FOUND)
  set(CMAKE_REQUIRED_LIBRARIES ${ABACUS_DEEPMD_TARGET})
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
  set(DeePMD_INPUT_NLIST_COMM_ARGS 0)
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

  # A communication constructor alone does not establish that the PyTorch
  # operator exchanges embeddings through MPI. Inspect the operator beside the
  # selected API library.
  get_target_property(_abacus_dp_api_library ${ABACUS_DEEPMD_TARGET} LOCATION)
  get_filename_component(_abacus_dp_library_dir "${_abacus_dp_api_library}"
                         DIRECTORY)
  unset(DeePMD_PT_OPERATOR CACHE)
  find_file(
    DeePMD_PT_OPERATOR libdeepmd_op_pt.so
    PATHS "${_abacus_dp_library_dir}"
    NO_DEFAULT_PATH)
  set(DeePMD_PT_OPERATOR_MPI FALSE)
  if(DeePMD_PT_OPERATOR AND CMAKE_READELF)
    execute_process(
      COMMAND "${CMAKE_READELF}" -d "${DeePMD_PT_OPERATOR}"
      RESULT_VARIABLE _abacus_dp_readelf_result
      OUTPUT_VARIABLE _abacus_dp_operator_dynamic)
    if(_abacus_dp_readelf_result EQUAL 0 AND _abacus_dp_operator_dynamic
                                             MATCHES "NEEDED[^\n]*libmpi\\.so")
      set(DeePMD_PT_OPERATOR_MPI TRUE)
    endif()
  endif()
  if(NOT DeePMD_PT_OPERATOR_MPI)
    set(DeePMD_INPUT_NLIST_COMM_ARGS 0)
  endif()

  if(ENABLE_MPI AND DEEPMD_BACKEND STREQUAL "Torch")
    if(NOT DeePMD_PT_OPERATOR_MPI)
      message(
        WARNING
          "Distributed PyTorch DeePMD inference is disabled: the operator's "
          "direct libmpi.so dependency could not be verified.")
    elseif(DeePMD_INPUT_NLIST_COMM_ARGS EQUAL 0)
      message(
        WARNING
          "Distributed PyTorch DeePMD inference is disabled: communication-capable "
          "InputNlist is unavailable.")
    endif()
  endif()
endif()

cmake_pop_check_state()
unset(_abacus_dp_check)
unset(_abacus_dp_nlist_type)
unset(_abacus_dp_nlist_header)
unset(_abacus_dp_api_library)
unset(_abacus_dp_library_dir)
unset(_abacus_dp_readelf_result)
unset(_abacus_dp_operator_dynamic)
mark_as_advanced(DeePMD_INCLUDE_DIR DeePMD_INCLUDE_C_DIR deepmd_cc deepmd_c
                 DeePMD_PT_OPERATOR)
