AddTest(
  TARGET hsolver_diago_elpa
  LIBS ${math_libs} ELPA::ELPA
  SOURCES diago_elpa_test.cpp
)

AddTest(
  TARGET hsolver_diago
  LIBS ${math_libs} ELPA::ELPA base orb
  SOURCES diago_test.cpp ../pdiag_double.cpp ../diag_scalapack_gvx.cpp 
          ../pdiag_common.cpp ../../src_parallel/parallel_global.cpp 
          ../../src_parallel/parallel_common.cpp ../../src_parallel/parallel_reduce.cpp 
)

install(FILES H-KPoints.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES H-GammaOnly.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES S-KPoints.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES S-GammaOnly.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES H-KPoints-large.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES H-GammaOnly-large.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES S-KPoints-large.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES S-GammaOnly-large.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES diago_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES diago_elpa_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(BASH bash)
add_test(NAME hsolver_diago_elpa_parallel
      COMMAND ${BASH} diago_elpa_parallel_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
add_test(NAME hsolver_diago_parallel
      COMMAND ${BASH} diago_parallel_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )    
