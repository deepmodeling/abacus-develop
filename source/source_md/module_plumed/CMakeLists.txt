add_library(
    module_plumed
    OBJECT
    plumed_interface.cpp
)

# The PLUMED interface is only compiled in when the build is configured with
# -DENABLE_PLUMED=ON
if(ENABLE_PLUMED)
  target_compile_definitions(module_plumed PRIVATE USE_PLUMED)
endif()

if(ENABLE_COVERAGE)
  add_coverage(module_plumed)
endif()
