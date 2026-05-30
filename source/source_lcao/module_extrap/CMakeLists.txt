list(APPEND objects
    wf_history_lcao.cpp
    wf_orthonormalize_lcao.cpp
)

add_library(
    lcao_extrap
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(lcao_extrap)
endif()
