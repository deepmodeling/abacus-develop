set(vdw_sources
    vdw_xcname.cpp
    vdwd2_parameters.cpp
    vdwd2.cpp
    vdwd3_data.cpp
    vdwd3_parameters.cpp
    vdwd3_evaluator.cpp
    vdwd3.cpp
    vdw.cpp
)

if(ENABLE_DFTD4)
  list(APPEND vdw_sources vdwd4.cpp)
endif()

add_library(vdw OBJECT ${vdw_sources})

# Density-grid nonlocal XC has a different dependency and data model from the
# atom-pair D2/D3/D4 correction. Keep it in a separate object target so the
# existing pairwise-vdW unit tests do not acquire PW/XC link dependencies.
set(nonlocal_xc_sources nonlocal_xc.cpp)

if(ENABLE_LIBVDWXC)
  list(APPEND nonlocal_xc_sources nonlocal_xc_libvdwxc.cpp)
endif()

add_library(nonlocal_xc OBJECT ${nonlocal_xc_sources})

if(ENABLE_DFTD4)
  target_link_libraries(vdw PUBLIC dftd4::dftd4)
endif()

if(ENABLE_COVERAGE)
  add_coverage(vdw)
endif()
