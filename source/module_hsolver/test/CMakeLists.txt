remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET HSolver_cg
  LIBS ${math_libs} base 
  SOURCES diago_cg_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp ../../module_hamilt/hamilt_pw.cpp 
          ../../module_psi/psi.cpp  ../../src_parallel/parallel_reduce.cpp
          ../../src_parallel/parallel_global.cpp ../../module_pw/unittest/test_tool.cpp
)
AddTest(
  TARGET HSolver_dav
  LIBS ${math_libs} base 
  SOURCES diago_david_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp ../../module_hamilt/hamilt_pw.cpp 
          ../../module_psi/psi.cpp  ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_global.cpp 
          ../../module_pw/unittest/test_tool.cpp
)
AddTest(
  TARGET HSolver_LCAO
  LIBS ${math_libs} ELPA::ELPA base
  SOURCES diago_lcao_test.cpp ../diago_elpa.cpp ../diago_blas.cpp ../../src_parallel/parallel_global.cpp 
          ../../src_parallel/parallel_common.cpp ../../src_parallel/parallel_reduce.cpp
)

install(FILES H-KPoints-Si2.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES H-GammaOnly-Si2.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES S-KPoints-Si2.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES S-GammaOnly-Si2.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES H-KPoints-Si64.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES H-GammaOnly-Si64.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES S-KPoints-Si64.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES S-GammaOnly-Si64.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES diago_cg_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES diago_david_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES diago_lcao_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(BASH bash)
add_test(NAME HSolver_cg_parallel
      COMMAND ${BASH} diago_cg_parallel_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
add_test(NAME HSolver_dav_parallel
      COMMAND ${BASH} diago_david_parallel_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    ) 
add_test(NAME HSolver_LCAO_parallel
      COMMAND ${BASH} diago_lcao_parallel_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )   