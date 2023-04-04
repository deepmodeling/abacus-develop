remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__MPI)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET io_rho_io
  LIBS ${math_libs} base device
  SOURCES rho_io_test.cpp ../read_rho.cpp ../write_rho.cpp ../output.cpp ../../module_cell/unitcell.cpp
  	../../module_cell/read_atoms.cpp ../../module_cell/read_cell_pseudopots.cpp ../../module_cell/atom_spec.cpp
	../../module_cell/atom_pseudo.cpp ../../module_cell/pseudo_nc.cpp ../../module_cell/read_pp.cpp
	../../module_cell/read_pp_upf201.cpp ../../module_cell/read_pp_upf100.cpp ../../module_cell/read_pp_vwr.cpp
	../../module_cell/read_pp_blps.cpp
)
