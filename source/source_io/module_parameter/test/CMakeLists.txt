# The availability parser and metadata validator depend only on the standard
# library. Keep this target outside ENABLE_MPI so serial-only builds exercise it.
AddTest(
  TARGET MODULE_IO_availability
  SOURCES availability_test.cpp
          ../availability.cpp
          ../availability_validator.cpp
)
