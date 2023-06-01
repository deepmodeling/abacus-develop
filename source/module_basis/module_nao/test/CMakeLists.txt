AddTest(
  TARGET numerical_radial
  SOURCES 
    numerical_radial_test.cpp
	../numerical_radial.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET atomic_radials
  SOURCES 
    atomic_radials_test.cpp
	../atomic_radials.cpp
	../radial_set.cpp
    ../numerical_radial.cpp
  LIBS ${math_libs} device base
)

install(DIRECTORY test_files DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
