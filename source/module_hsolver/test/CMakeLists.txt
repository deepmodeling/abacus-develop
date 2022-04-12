remove_definitions(-D__MPI)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET HSolver_cg
  LIBS ${math_libs} base 
  SOURCES diago_cg_test.cpp ../diagocg.cpp  ../iterdiagcon.cpp ../../module_hamilt/hamiltpw.cpp ../../module_psi/psi.cpp  ../../src_parallel/parallel_reduce.cpp
)

install(FILES data-H DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
