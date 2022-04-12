add_library(
    hsolver
    OBJECT
    diagocg.cpp
    diagodavid.cpp
    hsolverpw.cpp
    iterdiagcon.cpp
    hsolverlcao.cpp
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
