target_sources(
    md
    PRIVATE
    snap_calculator.cpp
    snap_descriptor.cpp
    snap_model.cpp
    snap_clebsch_gordan.cpp
)

if(BUILD_TESTING)
  add_subdirectory(unittest)
endif()
