target_sources(md PRIVATE snap_calculator.cpp snap_descriptor.cpp
                          snap_derivative.cpp snap_model.cpp snap_cg.cpp)

if(BUILD_TESTING)
  add_subdirectory(unittest)
endif()
