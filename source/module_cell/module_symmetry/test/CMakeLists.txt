remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET symmetry_analysis
  LIBS base ${math_libs} device symmetry
  SOURCES symmetry_test.cpp symmetry_test_analysis.cpp
)
