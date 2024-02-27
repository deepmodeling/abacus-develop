remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)
remove_definitions(-DUSE_PAW)

AddTest(
  TARGET io_json_output_json
  LIBS ${math_libs} base device cell_info
  SOURCES para_json_test.cpp ../general_info.cpp ../init_info.cpp 
  ../../para_json.cpp ../abacusjson.cpp ../../input.cpp ../../output.cpp 
)