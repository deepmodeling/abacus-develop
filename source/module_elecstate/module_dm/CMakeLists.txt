if(ENABLE_LCAO)
  list(APPEND objects
      dmk_io.cpp
      density_matrix.cpp
  )
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()