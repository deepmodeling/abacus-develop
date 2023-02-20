remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET cell_read_pp
  LIBS ${math_libs} base device
  SOURCES read_pp_test.cpp ../read_pp.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp ../read_pp_vwr.cpp ../read_pp_blps.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp ../../module_io/output.cpp
)

AddTest(
  TARGET cell_pseudo_nc
  LIBS ${math_libs} base device
  SOURCES pseudo_nc_test.cpp ../pseudo_nc.cpp ../read_pp.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp ../read_pp_vwr.cpp ../read_pp_blps.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp ../../module_io/output.cpp

)

AddTest(
  TARGET cell_atom_pseudo
  LIBS ${math_libs} base device
  SOURCES atom_pseudo_test.cpp ../atom_pseudo.cpp ../pseudo_nc.cpp ../read_pp.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp ../read_pp_vwr.cpp ../read_pp_blps.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp ../../module_io/output.cpp
)

install(FILES bcast_atom_pseudo_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)

add_test(NAME cell_bcast_atom_pseudo_test
      COMMAND ${BASH} bcast_atom_pseudo_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET klist_test
  LIBS ${math_libs} base device
  SOURCES klist_test.cpp ../klist.cpp 
  ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_kpoints.cpp 
  ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp 
  ../../module_io/output.cpp ../module_symmetry/symm_other.cpp
  ../module_symmetry/symmetry_basic.cpp ../module_symmetry/symmetry.cpp
)
