remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

AddTest(
  TARGET klist_test
  LIBS ${math_libs} base symmetry device
  SOURCES klist_test.cpp ../klist.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp ../../src_io/output.cpp
)

AddTest(
   TARGET SRC_PW_UTs
   LIBS ${math_libs} base device
   SOURCES force_multi_device_test.cpp stress_multi_device_test.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp
)

install(FILES KPT KPT1 KPT2 KPT4   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
