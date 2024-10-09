remove_definitions(-D__MPI)

AddTest(
  TARGET Sto_Tool_UTs
  LIBS parameter ${math_libs} psi base device 
  SOURCES ../sto_tool.cpp test_sto_tool.cpp
)