add_library(
    hsolver
    OBJECT
    diagocg.cpp
    diagodavid.cpp
    hsolverpw.cpp
    iterdiagcon.cpp
    hsolverlcao.cpp
    diagosca.cpp
    diagoelpa.cpp
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
