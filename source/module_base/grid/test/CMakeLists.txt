remove_definitions(-D__MPI)

AddTest(
  TARGET test_delley
  SOURCES test_delley.cpp
  ../delley.cpp
  ../../ylm.cpp
  ../../timer.cpp
  ../../tool_quit.cpp
  ../../global_variable.cpp
  ../../global_file.cpp
  ../../global_function.cpp
)

