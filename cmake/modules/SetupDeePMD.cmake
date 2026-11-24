# ABACUS-specific API selection and distributed-inference capabilities.
# Library discovery and usage requirements come from DeePMDConfig.cmake.
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state(RESET)

# Do not reuse results from a different DeePMD installation.
foreach(_abacus_dp_check IN
        ITEMS DeePMD_CXX_API_USABLE DeePMD_HAS_COMM_NLIST_13
              DeePMD_HAS_COMM_NLIST_12 DeePMD_HAS_NLIST_MAPPING)
  unset(${_abacus_dp_check} CACHE)
endforeach()

if(TARGET DeePMD::deepmd_cc)
  # TODO: Remove when DeePMD keeps its build-only common.h PCH private.
  set_property(TARGET DeePMD::deepmd_cc PROPERTY INTERFACE_PRECOMPILE_HEADERS
                                                 "")
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

if(NOT ABACUS_DEEPMD_TARGET)
  message(
    FATAL_ERROR "DeePMD must provide a usable C++ API or the C API target.")
endif()

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
  if(_abacus_dp_readelf_result EQUAL 0 AND _abacus_dp_operator_dynamic MATCHES
                                           "NEEDED[^\n]*libmpi\\.so")
    set(DeePMD_PT_OPERATOR_MPI TRUE)
  endif()
endif()
if(NOT DeePMD_PT_OPERATOR_MPI)
  set(DeePMD_INPUT_NLIST_COMM_ARGS 0)
endif()

if(ENABLE_MPI AND DeePMD_PT_OPERATOR)
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

cmake_pop_check_state()
unset(_abacus_dp_check)
unset(_abacus_dp_nlist_type)
unset(_abacus_dp_nlist_header)
unset(_abacus_dp_api_library)
unset(_abacus_dp_library_dir)
unset(_abacus_dp_readelf_result)
unset(_abacus_dp_operator_dynamic)
mark_as_advanced(DeePMD_PT_OPERATOR)
