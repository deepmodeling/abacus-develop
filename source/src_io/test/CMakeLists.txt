remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET read_wfc_pw
  LIBS ${math_libs} base planewave
  SOURCES read_wfc_pw_test.cpp ../wf_io.cpp ../rwstream.cpp ../../src_parallel/parallel_reduce.cpp 
  	../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_common.cpp  
	../../input.cpp klist_mock.cpp
)

install(FILES WAVEFUNC1.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES WAVEFUNC2.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
