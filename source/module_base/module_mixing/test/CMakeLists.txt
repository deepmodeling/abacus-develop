remove_definitions(-D__MPI)
AddTest(
  TARGET base_broyden_mixing
  LIBS base device ${math_libs}
  SOURCES broyden_mixing_test.cpp
)
AddTest(
  TARGET base_pulay_mixing
  LIBS base device ${math_libs}
  SOURCES pulay_mixing_test.cpp
)
AddTest(
  TARGET base_plain_mixing
  LIBS ${math_libs}
  SOURCES plain_mixing_test.cpp ../mixing_data.cpp
)