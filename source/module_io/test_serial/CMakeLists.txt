remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__MPI)

add_library(
    io_input_serial
    OBJECT
    ../read_input_item_system.cpp
    ../read_input_item_elec_stru.cpp
    ../read_input_item_relax.cpp
    ../read_input_item_md.cpp
    ../read_input_item_ofdft.cpp
    ../read_input_item_sdft.cpp
    ../read_input_item_tddft.cpp
    ../read_input_item_deepks.cpp
    ../read_input_item_model.cpp
    ../read_input_item_postprocess.cpp  
    ../read_input_item_exx_dftu.cpp
    ../read_input_item_other.cpp   
    ../read_input_item_output.cpp
    ../read_input.cpp
    ../read_set_globalv.cpp
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
AddTest(
  TARGET io_read_input_serial
  LIBS parameter ${math_libs} io_input_serial 
  SOURCES read_input_test.cpp 
          ../../module_base/test/tool_quit_no_exit.cpp 
          ../../module_base/module_device/device.cpp
)

AddTest(
  TARGET io_read_item_serial
  LIBS parameter ${math_libs} base device io_input_serial 
  SOURCES read_input_item_test.cpp
)

AddTest(
  TARGET io_read_input_tool
  SOURCES read_input_tool_test.cpp
)

AddTest(
  TARGET io_rho_io
  LIBS parameter ${math_libs} base device cell_info 
  SOURCES rho_io_test.cpp ../read_cube.cpp ../write_cube.cpp ../output.cpp
)

AddTest(
  TARGET io_nscf_band
  LIBS parameter ${math_libs} base device
  SOURCES nscf_band_test.cpp ../nscf_band.cpp
)

AddTest(
  TARGET io_system_variable_test
  LIBS parameter ${math_libs} base device io_input_serial 
  SOURCES io_system_variable_test.cpp
)