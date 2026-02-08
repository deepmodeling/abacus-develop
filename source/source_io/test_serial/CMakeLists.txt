remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__MPI)

add_library(
    io_input_serial
    OBJECT
    ../module_parameter/read_input_item_system.cpp
    ../module_parameter/read_input_item_elec_stru.cpp
    ../module_parameter/read_input_item_relax.cpp
    ../module_parameter/read_input_item_md.cpp
    ../module_parameter/read_input_item_ofdft.cpp
    ../module_parameter/read_input_item_sdft.cpp
    ../module_parameter/read_input_item_tddft.cpp
    ../module_parameter/read_input_item_deepks.cpp
    ../module_parameter/read_input_item_model.cpp
    ../module_parameter/read_input_item_postprocess.cpp  
    ../module_parameter/read_input_item_exx_dftu.cpp
    ../module_parameter/read_input_item_other.cpp   
    ../module_parameter/read_input_item_output.cpp
    ../module_parameter/read_input.cpp
    ../module_parameter/read_set_globalv.cpp
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
AddTest(
  TARGET MODULE_IO_read_input_serial
  LIBS parameter ${math_libs} io_input_serial 
  SOURCES read_input_test.cpp 
          ../../source_base/test/tool_quit_no_exit.cpp 
          ../../source_base/module_device/device.cpp
          ../../source_base/module_device/output_device.cpp
)

AddTest(
  TARGET MODULE_IO_read_item_serial
  LIBS parameter ${math_libs} base device io_input_serial 
  SOURCES read_input_item_test.cpp
)

AddTest(
  TARGET MODULE_IO_read_input_tool
  SOURCES read_input_tool_test.cpp
)

AddTest(
  TARGET MODULE_IO_rho_io
  LIBS parameter ${math_libs} base device cell_info 
  SOURCES rho_io_test.cpp ../module_output/read_cube.cpp ../module_output/write_cube.cpp ../module_output/output.cpp
)

AddTest(
  TARGET MODULE_IO_write_bands
  LIBS parameter ${math_libs} base device
  SOURCES write_bands_test.cpp ../module_energy/write_bands.cpp
)

AddTest(
  TARGET MODULE_IO_system_variable_test
  LIBS parameter ${math_libs} base device io_input_serial 
  SOURCES io_system_variable_test.cpp
)
