remove_definitions(-D__MPI)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)
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


AddTest(
  TARGET hsolver_sum_band
  LIBS ${math_libs} base 
  SOURCES sum_band_test.cpp ../charge.cpp ../../src_parallel/parallel_reduce.cpp ../../src_parallel/ft.cpp ../pw_basis.cpp ../pw_complement.cpp ../klist.cpp ../../src_parallel/parallel_kpoints.cpp 
)

install(FILES STRU DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES H.pz-vbc.UPF DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES KPT DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES SPIN1_CHG DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES WAVEFUNC1.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
