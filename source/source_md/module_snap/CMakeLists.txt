target_sources(
    md
    PRIVATE
    calc.cpp
    desc.cpp
    deriv.cpp
    model.cpp
    cg.cpp
)

if(BUILD_TESTING)
  add_subdirectory(unittests)
endif()
