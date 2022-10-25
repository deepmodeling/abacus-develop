list(APPEND objects
    parallel_common.cpp
    parallel_global.cpp
    parallel_grid.cpp
    parallel_kpoints.cpp
    parallel_reduce.cpp
)
if(ENABLE_LCAO)
  list(APPEND objects
    subgrid_oper.cpp
  )
endif()

add_library(
    parallel
    OBJECT
    ${objects}
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
