if (ENABLE_LIBRI)
  if(ENABLE_LCAO)
    add_library(
            rpa
            OBJECT
            rpa.cpp
    )

    if(ENABLE_COVERAGE)
      add_coverage(rpa)
    endif()
  endif()
endif ()
