add_library(
    lcao_extrap
    OBJECT
    wf_history_lcao.cpp
    wf_orthonormalize_lcao.cpp
)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

if(ENABLE_COVERAGE)
  add_coverage(lcao_extrap)
endif()
