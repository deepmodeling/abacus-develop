add_definitions(-D__NORMAL)
if (USE_CUDA)
AddTest(
  TARGET pw_test_gpu
  LIBS parameter ${math_libs} planewave device FFTW3::FFTW3_FLOAT
  SOURCES ../../../module_base/matrix.cpp ../../../module_base/complexmatrix.cpp ../../../module_base/matrix3.cpp ../../../module_base/tool_quit.cpp
          ../../../module_base/mymath.cpp ../../../module_base/timer.cpp ../../../module_base/memory.cpp ../../../module_base/blas_connector.cpp 
          ../../../module_base/libm/branred.cpp ../../../module_base/libm/sincos.cpp 
          ../../../module_base/module_device/memory_op.cpp
          depend_mock.cpp pw_test.cpp pw_basis_C2C.cpp pw_basis_C2R.cpp pw_basis_k_C2C.cpp 
          test_tool.cpp 
)
endif()

