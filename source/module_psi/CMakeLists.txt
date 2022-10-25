add_library(
    psi
    OBJECT
    psi.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(psi)
endif()
