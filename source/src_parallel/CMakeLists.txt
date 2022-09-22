list(APPEND objects
    parallel_common.cpp
    parallel_global.cpp
    parallel_grid.cpp
    parallel_kpoints.cpp
    parallel_reduce.cpp
    subgrid_oper.cpp
)

add_library(
    parallel
    OBJECT
    ${objects}
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
