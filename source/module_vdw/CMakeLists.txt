add_library(
    base
    OBJECT
    xxx.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(base)
endif()
