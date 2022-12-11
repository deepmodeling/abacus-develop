list(APPEND objects
      dftu.cpp
      dftu_relax.cpp
      dftu_yukawa.cpp
)

add_library(
    dftu
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(dftu)
endif()
