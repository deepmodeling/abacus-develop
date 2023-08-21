if(ENABLE_LCAO)
  list(APPEND objects
      density_matrix.cpp
  )
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()