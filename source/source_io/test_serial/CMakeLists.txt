abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__MPI)

add_library(
  io_input_serial OBJECT
  ../module_parameter/read_inp_sys.cpp
  ../module_parameter/read_inp_estruc.cpp
  ../module_parameter/read_inp_relax.cpp
  ../module_parameter/read_input_item_md.cpp
  ../module_parameter/read_inp_ofdft.cpp
  ../module_parameter/read_input_item_sdft.cpp
  ../module_parameter/read_inp_tddft.cpp
  ../module_parameter/read_inp_bse.cpp
  ../module_parameter/read_inp_dfpt.cpp
  ../module_parameter/read_inp_deepks.cpp
  ../module_parameter/read_inp_model.cpp
  ../module_parameter/read_inp_sccs.cpp
  ../module_parameter/sccs_inp_check.cpp
  ../module_parameter/read_inp_postproc.cpp
  ../module_parameter/read_inp_exx_dftu.cpp
  ../module_parameter/read_inp_other.cpp
  ../module_parameter/read_inp_out.cpp
  ../module_parameter/read_input.cpp
  ../module_parameter/read_set_globalv.cpp
  ../module_parameter/availability.cpp
  ../module_parameter/availability_validator.cpp
  ../module_parameter/input_domain.cpp
  # read_set_globalv.cpp references stru_frame_io::extract_last_frame, so
  # every target linking io_input_serial needs this translation unit.
  ../../source_cell/stru_frame_io.cpp)

abacus_copy_test_files(support)

addtest(
  TARGET MODULE_IO_read_input_serial
  LIBS parameter io_input_serial
  SOURCES read_input_test.cpp ../../source_base/test/tool_quit_no_exit.cpp
          ../../source_base/module_device/device.cpp
          ../../source_base/module_device/output_device.cpp)

# The serial-reader targets stay in this directory, which disables __MPI for
# their translation units, while their sources live in the sibling
# `unittests/` directory required by the test layout rule.
addtest(
  TARGET MODULE_IO_read_item_serial
  LIBS parameter base device io_input_serial
  SOURCES ../unittests/test_read_input_item.cpp test_read_inp_sccs.cpp)

addtest(
  TARGET MODULE_IO_input_domain_reader
  LIBS parameter base device io_input_serial
  SOURCES ../unittests/test_input_domain_reader.cpp ../input_help.cpp)

addtest(TARGET MODULE_IO_read_input_tool SOURCES read_input_tool_test.cpp)

addtest(
  TARGET MODULE_IO_rho_io
  LIBS parameter base device cell_info
  SOURCES ../unittests/test_rho_io.cpp ../module_output/read_cube.cpp
          ../module_output/write_cube.cpp)

addtest(
  TARGET MODULE_IO_read_cube
  LIBS parameter base device cell_info
  SOURCES test_read_cube.cpp ../module_output/read_cube.cpp)

addtest(
  TARGET MODULE_IO_system_variable_test
  LIBS parameter base device io_input_serial
  SOURCES ../unittests/test_io_sys_var.cpp)
