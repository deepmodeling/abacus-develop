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
  add_coverage(nonlocal_xc)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
  if(ENABLE_LIBVDWXC AND NOT ENABLE_MPI AND NOT USE_CUDA AND NOT USE_ROCM)
    # The serial adapter test is registered here because test/ is MPI-only.
    AddTest(
      TARGET MODULE_HAMILT_vdwTest_LIBVDWXC
      LIBS xc_ planewave base parameter device nonlocal_xc
      SOURCES test/test_nonlocal_xc_libvdwxc.cpp
    )
    target_link_libraries(MODULE_HAMILT_vdwTest_LIBVDWXC PRIVATE Libxc::xc libvdwxc::vdwxc)
    set_tests_properties(MODULE_HAMILT_vdwTest_LIBVDWXC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
  endif()
endif()
