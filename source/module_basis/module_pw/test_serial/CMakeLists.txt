remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)
AddTest(
  TARGET basis_pw_serial
  LIBS ${math_libs}
  SOURCES pw_basis_test.cpp ../pw_basis.cpp ../pw_init.cpp
  		../pw_distributeg.cpp
  		../pw_distributer.cpp
		../pw_distributeg_method1.cpp
		../pw_distributeg_method2.cpp
  		../fft.cpp
		../../../module_base/blas_connector.h
		../../../module_base/matrix.cpp
		../../../module_base/matrix3.cpp
		../../../module_base/timer.cpp
		../../../module_base/mymath.cpp
		../../../module_base/tool_quit.cpp
		../../../module_base/memory.cpp
		../../../module_base/global_variable.cpp
		../../../module_base/global_function.cpp
		../../../module_base/global_file.cpp
)
