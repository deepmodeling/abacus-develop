remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__UT_USE_CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__UT_USE_ROCM)
remove_definitions(-D__DEEPKS)

add_library(
    planewave_serial
    OBJECT
    ../fft.cpp
    ../pw_basis.cpp
    ../pw_basis_k.cpp
    ../pw_basis_sup.cpp
    ../pw_distributeg.cpp
    ../pw_distributeg_method1.cpp
    ../pw_distributeg_method2.cpp
    ../pw_distributer.cpp
    ../pw_init.cpp
    ../pw_transform.cpp
    ../pw_transform_k.cpp
)

AddTest(
  TARGET basis_pw_serial
  LIBS ${math_libs} planewave_serial device base
  SOURCES pw_basis_test.cpp
)

AddTest(
  TARGET basis_pw_k_serial
  LIBS ${math_libs} planewave_serial device base
  SOURCES pw_basis_k_test.cpp
)
