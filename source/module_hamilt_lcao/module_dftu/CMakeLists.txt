if(ENABLE_LCAO)
  list(APPEND objects
        dftu.cpp
        dftu_force.cpp
        dftu_yukawa.cpp
        dftu_folding.cpp
        dftu_io.cpp
        dftu_tools.cpp
        dftu_occup.cpp
        dftu_hamilt.cpp
  )

  add_library(
      dftu
      OBJECT
      ${objects}
  )

  if(ENABLE_COVERAGE)
    add_coverage(dftu)
  endif()
endif()