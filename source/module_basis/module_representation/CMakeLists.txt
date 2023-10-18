list(APPEND objects
    representation.cpp
    repin.cpp
    repin_pao.cpp
    repin_nao.cpp
    repout.cpp
    repout_pw.cpp
)

add_library(
    representation
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(representation)
endif()

#if(BUILD_TESTING)
#  add_subdirectory(test)
#endif()
