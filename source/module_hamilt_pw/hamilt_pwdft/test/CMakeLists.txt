remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

AddTest(
  TARGET pwdft_soc
  LIBS parameter  ${math_libs}
  SOURCES soc_test.cpp ../soc.cpp
  	../../../module_base/global_variable.cpp
  	../../../module_base/global_function.cpp
	../../../module_base/global_file.cpp
	../../../module_base/tool_quit.cpp
	../../../module_base/global_variable.cpp
	../../../module_base/global_file.cpp
	../../../module_base/memory.cpp
	../../../module_base/timer.cpp
	../../../module_base/blas_connector.cpp
	../../../module_base/parallel_global.cpp
	../../../module_base/parallel_comm.cpp
	../../../module_base/parallel_common.cpp
	../../../module_base/parallel_reduce.cpp
)

AddTest(
	TARGET radial_proj_test
	LIBS parameter  base device ${math_libs}
	SOURCES radial_proj_test.cpp ../radial_proj.cpp
)