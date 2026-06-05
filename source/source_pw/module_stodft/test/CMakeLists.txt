remove_definitions(-D__MPI)

AddTest(
  TARGET MODULE_PW_Sto_Tool_UTs
  LIBS parameter ${math_libs} psi base device 
  SOURCES ../sto_tool.cpp test_sto_tool.cpp
)

AddTest(
  TARGET MODULE_PW_Sto_Hamilt_UTs
  LIBS parameter ${math_libs} psi base device planewave_serial 
  SOURCES ../hamilt_sdft_pw.cpp test_hamilt_sto.cpp ../../../source_hamilt/operator.cpp
)