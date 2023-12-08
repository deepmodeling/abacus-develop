remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET relax_new_line_search
  LIBS formatter
  SOURCES line_search_test.cpp ../line_search.cpp ../../../module_base/global_variable.cpp ../../../module_base/global_file.cpp ../../../module_base/global_function.cpp ../../../module_base/memory.cpp ../../../module_base/timer.cpp ../../../module_base/tool_quit.cpp
)

AddTest(
  TARGET relax_new_relax
  SOURCES relax_test.cpp ../relax.cpp ../line_search.cpp ../../../module_base/tool_quit.cpp ../../../module_base/global_variable.cpp ../../../module_base/global_file.cpp ../../../module_base/memory.cpp ../../../module_base/timer.cpp
    ../../../module_base/matrix3.cpp ../../../module_base/intarray.cpp ../../../module_base/tool_title.cpp ../../../module_io/input.cpp
    ../../../module_base/global_function.cpp ../../../module_base/complexmatrix.cpp ../../../module_base/matrix.cpp
    ../../../module_base/complexarray.cpp ../../../module_base/tool_quit.cpp ../../../module_base/realarray.cpp
  LIBS ${math_libs} formatter
)