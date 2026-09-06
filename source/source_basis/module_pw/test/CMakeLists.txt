abacus_add_local_feature_definitions(__NORMAL)
AddTest(
  TARGET MODULE_PW_pw_test
  LIBS parameter planewave device
  SOURCES ../../../source_base/matrix.cpp ../../../source_base/complexmatrix.cpp ../../../source_base/matrix3.cpp ../../../source_base/tool_quit.cpp
          ../../../source_base/mymath.cpp ../../../source_base/timer.cpp ../../../source_base/memory_recorder.cpp
          ../../../source_base/module_external/blas_connector_base.cpp ../../../source_base/module_external/blas_connector_vector.cpp ../../../source_base/module_external/blas_connector_matrix.cpp   
          ../../../source_base/libm/branred.cpp ../../../source_base/libm/sincos.cpp 
          ../../../source_base/module_device/memory_op.cpp
          depend_mock.cpp pw_test.cpp test1_1_1.cpp test1_1_2.cpp test1_2.cpp test1_3.cpp test1_4.cpp  test1_5.cpp
          test2_1_1.cpp test2_1_2.cpp test2_2.cpp test2_3.cpp
          test3_1.cpp test3_2.cpp test3_3.cpp test3_3_2.cpp
          test4_1.cpp test4_2.cpp test4_3.cpp test4_4.cpp  test4_5.cpp
          test5_1_1.cpp test5_1_2.cpp test5_2_1.cpp test5_2_2.cpp test5_3_1.cpp test5_4_1.cpp test5_4_2.cpp
          test6_1_1.cpp test6_1_2.cpp test6_2_1.cpp test6_2_2.cpp test6_3_1.cpp test6_4_1.cpp test6_4_2.cpp
          test7_1.cpp test6_2_1.cpp test7_3_1.cpp test7_3_2.cpp
          test8_1.cpp test8_2_1.cpp test8_3_1.cpp test8_3_2.cpp
          test_tool.cpp test_big.cpp test_other.cpp test_sup.cpp
)

add_test(NAME MODULE_PW_pw_test_parallel
      COMMAND mpirun -np 3 ./MODULE_PW_pw_test; mpirun -np 4 ./MODULE_PW_pw_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
