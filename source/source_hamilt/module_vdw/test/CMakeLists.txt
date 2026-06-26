remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/c6.txt ${CMAKE_CURRENT_BINARY_DIR}/c6.txt COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/r0.txt ${CMAKE_CURRENT_BINARY_DIR}/r0.txt COPYONLY)

AddTest(
  TARGET MODULE_HAMILT_vdwTest
  LIBS parameter ${math_libs} base device vdw 
  SOURCES vdw_test.cpp 
)

if(ENABLE_DFTD4)
  target_compile_definitions(MODULE_HAMILT_vdwTest PRIVATE __DFTD4)
  target_link_libraries(MODULE_HAMILT_vdwTest dftd4::dftd4)
endif()
