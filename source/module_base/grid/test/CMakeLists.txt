remove_definitions(-D__MPI)

AddTest(
  TARGET test_delley
  SOURCES test_delley.cpp
  ../delley.cpp
  ../../ylm.cpp
)

AddTest(
  TARGET test_radial
  SOURCES test_radial.cpp
  ../radial.cpp
)
