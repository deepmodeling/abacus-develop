# The availability parser and metadata validator depend only on the standard
# library. Keep this target outside ENABLE_MPI so serial-only builds exercise it.
addtest(TARGET MODULE_IO_availability
        SOURCES availability_test.cpp ../availability.cpp
                ../availability_validator.cpp)

addtest(TARGET MODULE_IO_input_domain
        SOURCES test_input_domain.cpp ../input_domain.cpp ../availability.cpp
                ../availability_validator.cpp)
