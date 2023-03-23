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
  TARGET io_dos_test
  LIBS ${math_libs} base device
  SOURCES dos_test.cpp ../dos.cpp
)

AddTest(
  TARGET io_write_dos_pw
  LIBS ${math_libs} base device symmetry
  SOURCES write_dos_pw_test.cpp ../dos.cpp ../write_dos_pw.cpp ../output.cpp ../../module_cell/parallel_kpoints.cpp ../../module_cell/klist.cpp
)

AddTest(
  TARGET io_print_info
  LIBS ${math_libs} base device
  SOURCES print_info_test.cpp ../print_info.cpp ../input.cpp ../output.cpp ../../module_cell/unitcell.cpp
  	../../module_cell/read_atoms.cpp ../../module_cell/read_cell_pseudopots.cpp ../../module_cell/atom_spec.cpp
	../../module_cell/atom_pseudo.cpp ../../module_cell/pseudo_nc.cpp ../../module_cell/read_pp.cpp
	../../module_cell/read_pp_upf201.cpp ../../module_cell/read_pp_upf100.cpp ../../module_cell/read_pp_vwr.cpp
	../../module_cell/read_pp_blps.cpp ../../module_cell/klist.cpp ../../module_cell/parallel_kpoints.cpp 
	../../module_cell/module_symmetry/symm_other.cpp ../../module_cell/module_symmetry/symmetry_basic.cpp
	../../module_cell/module_symmetry/symmetry.cpp
)

 AddTest(
  TARGET io_write_input_test
  LIBS ${math_libs} base device
  SOURCES write_input_test.cpp ../write_input.cpp ../input.cpp
)
