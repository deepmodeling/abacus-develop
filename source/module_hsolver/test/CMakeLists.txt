remove_definitions(-D__MPI)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET HSolver_cg
  LIBS ${math_libs} base 
  SOURCES diago_cg_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp ../../module_hamilt/hamilt_pw.cpp ../../module_psi/psi.cpp  ../../src_parallel/parallel_reduce.cpp
)
AddTest(
  TARGET HSolver_dav
  LIBS ${math_libs} base 
  SOURCES diago_david_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp ../../module_hamilt/hamilt_pw.cpp ../../module_psi/psi.cpp  ../../src_parallel/parallel_reduce.cpp
)

install(FILES data-H DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
