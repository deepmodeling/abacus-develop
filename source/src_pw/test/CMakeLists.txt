remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

AddTest(
  TARGET klist_test
  LIBS ${math_libs} base symmetry
  SOURCES klist_test.cpp ../klist.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp ../../src_io/output.cpp
)

install(FILES KPT KPT1 KPT2 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
