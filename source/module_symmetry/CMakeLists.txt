add_library(
    symmetry
    OBJECT
    symm_other.cpp
    symmetry_basic.cpp
    symmetry.cpp
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
