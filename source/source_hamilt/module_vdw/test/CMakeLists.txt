abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

install(FILES c6.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES r0.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_HAMILT_vdwTest
  LIBS parameter base device vdw
  SOURCES vdw_test.cpp 
)

if(ENABLE_DFTD4)
  target_compile_definitions(MODULE_HAMILT_vdwTest PRIVATE __DFTD4)
  target_link_libraries(MODULE_HAMILT_vdwTest PRIVATE dftd4::dftd4)
endif()
