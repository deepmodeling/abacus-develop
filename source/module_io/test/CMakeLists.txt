remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET io_input_test
  LIBS ${math_libs} base device
  SOURCES input_test.cpp ../input.cpp
)

AddTest(
  TARGET io_input_test_para
  LIBS ${math_libs} base device
  SOURCES input_test_para.cpp ../input.cpp
)

install(FILES input_test_para.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME io_input_test_para4
      COMMAND ${BASH} input_test_para.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET io_input_conv
  LIBS ${math_libs} base device
  SOURCES input_conv_test.cpp ../input.cpp ../input_conv.cpp
)

AddTest(
  TARGET io_winput_test
  LIBS ${math_libs} base device
  SOURCES winput_test.cpp ../winput.cpp
)

AddTest(
  TARGET io_output_test
  LIBS ${math_libs} base device
  SOURCES output_test.cpp ../output.cpp
)

AddTest(
  TARGET binstream_test
  SOURCES binstream_test.cpp ../binstream.cpp
)

AddTest(
  TARGET io_write_istate_info_test
  LIBS ${math_libs} base device symmetry
  SOURCES write_istate_info_test.cpp ../write_istate_info.cpp ../output.cpp ../../module_cell/parallel_kpoints.cpp ../../module_cell/klist.cpp
)

AddTest(
  TARGET io_cal_dos
  LIBS ${math_libs} base device
  SOURCES cal_dos_test.cpp ../cal_dos.cpp
)

AddTest(
  TARGET io_write_dos_pw
  LIBS ${math_libs} base device symmetry
  SOURCES write_dos_pw_test.cpp ../cal_dos.cpp ../write_dos_pw.cpp ../output.cpp ../../module_cell/parallel_kpoints.cpp ../../module_cell/klist.cpp
)

AddTest(
  TARGET io_print_info
  LIBS ${math_libs} base device symmetry cell_info
  SOURCES print_info_test.cpp ../print_info.cpp ../input.cpp ../output.cpp ../../module_cell/klist.cpp ../../module_cell/parallel_kpoints.cpp
)

AddTest(
  TARGET io_write_input_test
  LIBS ${math_libs} base device
  SOURCES write_input_test.cpp ../write_input.cpp ../input.cpp
)

AddTest(
  TARGET io_single_R_test
  LIBS ${math_libs}
  SOURCES single_R_io_test.cpp ../single_R_io.cpp
  	../../module_base/global_variable.cpp
	../../module_base/parallel_reduce.cpp
	../../module_base/parallel_common.cpp
	../../module_base/parallel_global.cpp
)

AddTest(
  TARGET io_write_wfc_nao
  LIBS ${math_libs} base device
  SOURCES write_wfc_nao_test.cpp ../write_wfc_nao.cpp
)

AddTest(
  TARGET io_write_orb_info
  LIBS ${math_libs} base device cell_info
  SOURCES write_orb_info_test.cpp ../write_orb_info.cpp ../output.cpp
)

AddTest(
  TARGET io_parse_args
  SOURCES parse_args_test.cpp
)

AddTest(
  TARGET io_bessel_basis_test
  LIBS ${math_libs} base device
  SOURCES bessel_basis_test.cpp ../bessel_basis.cpp
)

AddTest(
  TARGET io_output_log_test
  LIBS base ${math_libs} device
  SOURCES ../output_log.cpp outputlog_test.cpp
)

AddTest(
  TARGET io_sparse_matrix_test
  LIBS base ${math_libs} device
  SOURCES sparse_matrix_test.cpp ../sparse_matrix.cpp
)

AddTest(
  TARGET io_output_hcontainer_test
  LIBS base ${math_libs} device hcontainer
  SOURCES output_hcontainer_test.cpp ../output_hcontainer.cpp ../../module_basis/module_ao/parallel_2d.cpp ../../module_basis/module_ao/parallel_orbitals.cpp ../sparse_matrix.cpp tmp_mocks.cpp
)