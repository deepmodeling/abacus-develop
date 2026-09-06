abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__MPI)

add_library(
    io_input_serial
    OBJECT
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
    ../module_parameter/read_inp_postproc.cpp  
    ../module_parameter/read_inp_exx_dftu.cpp
    ../module_parameter/read_inp_other.cpp   
    ../module_parameter/read_inp_out.cpp
    ../module_parameter/read_input.cpp
    ../module_parameter/read_set_globalv.cpp
    ../module_parameter/availability.cpp
    ../module_parameter/availability_validator.cpp
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
AddTest(
  TARGET MODULE_IO_read_input_serial
  LIBS parameter io_input_serial
  SOURCES read_input_test.cpp 
          ../../source_base/test/tool_quit_no_exit.cpp 
          ../../source_base/module_device/device.cpp
          ../../source_base/module_device/output_device.cpp
)

AddTest(
  TARGET MODULE_IO_read_item_serial
  LIBS parameter base device io_input_serial
  SOURCES read_input_item_test.cpp
)

AddTest(
  TARGET MODULE_IO_read_input_tool
  SOURCES read_input_tool_test.cpp
)

AddTest(
  TARGET MODULE_IO_rho_io
  LIBS parameter base device cell_info
  SOURCES rho_io_test.cpp ../module_output/read_cube.cpp ../module_output/write_cube.cpp
)

AddTest(
  TARGET MODULE_IO_write_bands
  LIBS parameter base device symmetry
  SOURCES write_bands_test.cpp ../module_output/band_parallel_output.cpp ../module_energy/write_bands.cpp
  ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_IO_system_variable_test
  LIBS parameter base device io_input_serial
  SOURCES io_sys_var_test.cpp
)
