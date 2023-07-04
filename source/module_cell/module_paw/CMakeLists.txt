add_library(
    paw
    OBJECT
    paw_element.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(paw)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()