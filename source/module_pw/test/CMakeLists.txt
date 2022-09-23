add_definitions(-D__NORMAL)
AddTest(
  TARGET pw_test
  LIBS ${math_libs} planewave
  SOURCES ../../module_base/matrix.cpp ../../module_base/matrix3.cpp ../../module_base/tool_quit.cpp 
          ../../module_base/mymath3.cpp ../../module_base/timer.cpp ../../module_base/global_variable.cpp 
          ../../src_parallel/parallel_global.cpp
          pw_test.cpp test1-1-1.cpp test1-1-2.cpp test1-2.cpp test1-3.cpp test1-4.cpp  test1-5.cpp
          test2-1-1.cpp test2-1-2.cpp test2-2.cpp test2-3.cpp 
          test3-1.cpp test3-2.cpp test3-3.cpp 
          test4-1.cpp test4-2.cpp test4-3.cpp test4-4.cpp  test4-5.cpp
          test_tool.cpp
)

add_test(NAME pw_test_parallel
      COMMAND mpirun -np 2 ./pw_test; mpirun -np 3./pw_test
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
