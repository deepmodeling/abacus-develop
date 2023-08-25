list(APPEND objects
    basic_funcs.cpp
    lambda_loop.cpp
)

add_library(
    deltaspin
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(deltaspin)
endif()

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
