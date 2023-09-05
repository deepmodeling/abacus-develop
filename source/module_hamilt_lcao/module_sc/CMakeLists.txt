if(ENABLE_LCAO)
  list(APPEND objects
      spin_constrain.cpp
      sc_parse_json.cpp
      cal_weight_func.cpp
  )

  add_library(
      sc
      OBJECT
      ${objects}
  )

  if(ENABLE_COVERAGE)
    add_coverage(sc)
  endif()

  if(BUILD_TESTING)
      add_subdirectory(test)
  endif()
endif()

