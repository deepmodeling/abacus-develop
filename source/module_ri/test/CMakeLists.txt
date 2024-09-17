remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET dm_mixing_test
  LIBS parameter base ${math_libs} device
  SOURCES dm_mixing_test.cpp ../Mix_DMk_2D.cpp ../Mix_Matrix.cpp
)
AddTest(
  TARGET ri_cv_io_test
  LIBS base ${math_libs} device parameter
  SOURCES ri_cv_io_test.cpp
)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})