remove_definitions(-D__LCAO)
remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET MODULE_CELL_SYMMETRY_analysis
  LIBS parameter base ${math_libs} device symmetry 
  SOURCES symmetry_test.cpp symmetry_test_analysis.cpp
)
AddTest(
  TARGET MODULE_CELL_SYMMETRY_symtrz
  LIBS parameter base ${math_libs} device symmetry 
  SOURCES symmetry_test.cpp symmetry_test_symtrz.cpp
)