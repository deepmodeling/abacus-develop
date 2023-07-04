list(APPEND objects
    base_matrix.cpp
    atom_pair.cpp
    hcontainer.cpp
)

add_library(
    hcontainer
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hcontainer)
endif()

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()