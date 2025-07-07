remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

AddTest(
  TARGET MODULE_IO_JSON_OUTPUT_TEST
  LIBS parameter ${math_libs} base device cell_info 
  SOURCES para_json_test.cpp ../general_info.cpp ../init_info.cpp ../readin_info.cpp 
  ../../para_json.cpp ../abacusjson.cpp ../../output.cpp 
)