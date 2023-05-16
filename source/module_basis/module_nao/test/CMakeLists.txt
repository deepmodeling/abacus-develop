AddTest(
  TARGET numerical_radial
  SOURCES 
    numerical_radial_test.cpp
	../numerical_radial.cpp
  LIBS ${math_libs} device base
)

