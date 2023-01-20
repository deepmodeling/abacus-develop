remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET relax_new_line_search
  SOURCES line_search_test.cpp ../line_search.cpp ../../../module_base/tool_quit.cpp ../../../module_base/global_variable.cpp ../../../module_base/global_file.cpp ../../../module_base/memory.cpp ../../../module_base/timer.cpp
)
