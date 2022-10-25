list(APPEND objects
    elecstate.cpp
    elecstate_pw.cpp
    elecstate_pw_sdft.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      elecstate_lcao.cpp
      elecstate_lcao_tddft.cpp
      dm2d_to_grid.cpp
  )
endif()

add_library(
    elecstate
    OBJECT
    ${objects}
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
