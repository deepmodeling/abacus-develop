remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET HSolver_cg
  LIBS ${math_libs} base psi
  SOURCES diago_cg_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp  
          ../../src_parallel/parallel_reduce.cpp
          ../../src_parallel/parallel_global.cpp ../../module_pw/test/test_tool.cpp
          ../../module_hamilt/operator.cpp
          ../../module_hamilt/ks_pw/operator_pw.cpp
          ../src/math_kernel.cpp
          ../src/dngvd_op.cpp
)
AddTest(
  TARGET HSolver_dav
  LIBS ${math_libs} base psi
  SOURCES diago_david_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp 
          ../../src_parallel/parallel_reduce.cpp ../../src_parallel/parallel_global.cpp  ../../src_parallel/parallel_common.cpp
          ../../module_pw/test/test_tool.cpp
          ../../module_hamilt/operator.cpp
          ../../module_hamilt/ks_pw/operator_pw.cpp
          ../src/math_kernel.cpp
          ../src/dngvd_op.cpp
)

if(ENABLE_LCAO)
  if(USE_ELPA)
  AddTest(
    TARGET HSolver_LCAO
    LIBS ${math_libs} ELPA::ELPA base genelpa psi
    SOURCES diago_lcao_test.cpp ../diago_elpa.cpp ../diago_blas.cpp ../../src_parallel/parallel_global.cpp 
            ../../src_parallel/parallel_common.cpp ../../src_parallel/parallel_reduce.cpp
  )
  else()
    AddTest(
      TARGET HSolver_LCAO
      LIBS ${math_libs} base psi
      SOURCES diago_lcao_test.cpp ../diago_blas.cpp ../../src_parallel/parallel_global.cpp
      ../../src_parallel/parallel_common.cpp ../../src_parallel/parallel_reduce.cpp
    )
  endif()
endif()

if (USE_CUDA)
AddTest(
  TARGET Hsolver_UTs
  LIBS ${math_libs} base cublas
  SOURCES math_kernel_test.cpp ../../src_parallel/parallel_global.cpp 
          ../../src_parallel/parallel_common.cpp  ../../src_parallel/parallel_reduce.cpp
          ../../module_psi/src/memory_psi.cpp
          ../../module_hsolver/src/math_kernel.cpp
          ../../module_psi/src/cuda/memory.cu
          ../src/cuda/math_kernel.cu
)
AddTest(
  TARGET Dngvd_UTs
  LIBS ${math_libs} base cublas cusolver
  SOURCES math_dngvd_test.cpp 
          ../../src_parallel/parallel_global.cpp 
          ../../src_parallel/parallel_common.cpp  ../../src_parallel/parallel_reduce.cpp
          ../../module_psi/src/memory_psi.cpp
          ../../module_psi/src/cuda/memory.cu
          ../../module_hsolver/src/math_kernel.cpp
          ../src/cuda/math_kernel.cu
          ../../module_hsolver/src/dngvd_op.cpp
          ../src/cuda/dngvd_op.cu

)
else()
AddTest(
  TARGET Hsolver_UTs
  LIBS ${math_libs} base
  SOURCES math_kernel_test.cpp ../../src_parallel/parallel_global.cpp 
          ../../src_parallel/parallel_common.cpp  ../../src_parallel/parallel_reduce.cpp
          ../../module_psi/src/memory_psi.cpp
          ../src/math_kernel.cpp
)
endif()

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
if(ENABLE_LCAO)
  add_test(NAME HSolver_LCAO_parallel
        COMMAND ${BASH} diago_lcao_parallel_test.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
endif()   