remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

AddTest(
  TARGET cell_read_pp
  LIBS ${math_libs} base device
  SOURCES read_pp_test.cpp ../read_pp.cpp ../read_pp_upf201.cpp ../read_pp_upf100.cpp ../read_pp_vwr.cpp ../read_pp_blps.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp ../../module_io/output.cpp
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
