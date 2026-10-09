add_library(genelpa OBJECT
  elpa_new.cpp
  elpa_new_real.cpp
  elpa_new_complex.cpp
  elpa_runtime_check.cpp
  utils.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(genelpa)
endif()
