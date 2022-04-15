remove_definitions(-D__MPI)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET hsolver_david
  LIBS ${math_libs} base 
  SOURCES diago_david_test.cpp ../diago_david.cpp ../../src_parallel/parallel_reduce.cpp 
)

AddTest(
  TARGET hsolver_cg
  LIBS ${math_libs} base 
  SOURCES diago_cg_test.cpp ../diago_cg.cpp ../../src_parallel/parallel_reduce.cpp 
)

install(FILES data-H DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
