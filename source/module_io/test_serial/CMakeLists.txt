remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__MPI)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET io_rho_io
  LIBS ${math_libs} base device cell_simple
  SOURCES rho_io_test.cpp ../read_cube.cpp ../write_cube.cpp ../read_rho.cpp ../write_rho.cpp ../output.cpp
)

AddTest(
  TARGET io_dm_io
  LIBS ${math_libs} base device cell_simple
  SOURCES dm_io_test.cpp ../read_dm.cpp ../write_dm.cpp ../output.cpp
)

AddTest(
  TARGET io_nscf_band
  LIBS ${math_libs} base device
  SOURCES nscf_band_test.cpp ../nscf_band.cpp
)

AddTest(
  TARGET io_read_wfc_nao_test
  LIBS ${math_libs} base device
  SOURCES read_wfc_nao_test.cpp ../read_wfc_nao.cpp ../../module_psi/psi.cpp
)
