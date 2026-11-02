file(MAKE_DIRECTORY "${CASE_DIR}")
set(input_file "${CASE_DIR}/INPUT")
set(contents "INPUT_PARAMETERS\nbasis_type pw\nntype 1\n")
set(expected_exit 1)
set(expected_manual "input-main.html")

if(INPUT_CASE STREQUAL "missing_file")
  file(REMOVE "${input_file}")
  set(expected "Cannot open the INPUT file")
  set(expected_manual "")
elseif(INPUT_CASE STREQUAL "missing_header")
  set(contents "basis_type pw\nntype 1\n")
  set(expected "Missing INPUT_PARAMETERS header")
elseif(INPUT_CASE STREQUAL "unknown_parameter")
  string(APPEND contents "unknown_input_option 1\n")
  set(expected "Unknown INPUT parameter 'unknown_input_option'")
elseif(INPUT_CASE STREQUAL "duplicate_parameter")
  string(APPEND contents "nspin 1\nnspin 2\n")
  set(expected "The parameter nspin has been read twice")
  set(expected_manual "input-main.html#nspin")
elseif(INPUT_CASE STREQUAL "duplicate_string_eof")
  # Repeat a scalar string parameter with no value or final newline.
  string(APPEND contents "basis_type")
  set(expected "The parameter basis_type has been read twice")
  set(expected_manual "input-main.html#basis-type")
elseif(INPUT_CASE STREQUAL "invalid_value")
  string(APPEND contents "nspin 2abc\n")
  set(expected "Expected one Integer from {1, 2, 4}")
  set(expected_manual "input-main.html#nspin")
elseif(INPUT_CASE STREQUAL "invalid_dependency")
  string(APPEND contents "nspin 2\nnoncolin true\n")
  set(expected "Requires nspin==4")
  set(expected_manual "input-main.html#nspin")
elseif(INPUT_CASE STREQUAL "valid")
  string(APPEND contents "nspin 2\n")
  set(expected_exit 0)
else()
  message(FATAL_ERROR "Unknown INPUT test case: ${INPUT_CASE}")
endif()

if(NOT INPUT_CASE STREQUAL "missing_file")
  file(WRITE "${input_file}" "${contents}")
endif()

execute_process(
  COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          "${READER}" ${MPIEXEC_POSTFLAGS} "${input_file}"
  WORKING_DIRECTORY "${CASE_DIR}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error
  TIMEOUT 30)
string(APPEND output "${error}")

if(NOT "${result}" STREQUAL "${expected_exit}")
  message(
    FATAL_ERROR "Expected exit ${expected_exit}, got ${result}:\n${output}")
endif()
string(REGEX MATCHALL "INPUT reader process reached shutdown" exits "${output}")
list(LENGTH exits exit_count)
if(NOT exit_count EQUAL 4)
  message(
    FATAL_ERROR
      "Only ${exit_count} of 4 processes reached reader shutdown:\n${output}")
endif()
if(expected_exit EQUAL 1)
  foreach(diagnostic IN ITEMS "${input_file}" "${expected}"
                              "${expected_manual}")
    string(FIND "${output}" "${diagnostic}" found)
    if(found EQUAL -1)
      message(FATAL_ERROR "Missing diagnostic '${diagnostic}':\n${output}")
    endif()
  endforeach()
endif()
