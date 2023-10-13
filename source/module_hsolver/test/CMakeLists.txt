remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-DUSE_PAW)

AddTest(
  TARGET HSolver_bpcg
  LIBS ${math_libs} base psi device container
  SOURCES diago_bpcg_test.cpp ../diago_bpcg.cpp  ../diago_iter_assist.cpp  
          ../../module_basis/module_pw/test/test_tool.cpp
          ../../module_hamilt_general/operator.cpp
          ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
)
AddTest(
  TARGET HSolver_cg
  LIBS ${math_libs} base psi device
  SOURCES diago_cg_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp  
          ../../module_basis/module_pw/test/test_tool.cpp
          ../../module_hamilt_general/operator.cpp
          ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
)
AddTest(
  TARGET HSolver_cg_float
  LIBS ${math_libs} base psi device
  SOURCES diago_cg_float_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp  
          ../../module_basis/module_pw/test/test_tool.cpp
          ../../module_hamilt_general/operator.cpp
          ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
)
AddTest(
  TARGET HSolver_dav
  LIBS ${math_libs} base psi device
  SOURCES diago_david_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp 
          ../../module_basis/module_pw/test/test_tool.cpp
          ../../module_hamilt_general/operator.cpp
          ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
)
AddTest(
  TARGET HSolver_dav_float
  LIBS ${math_libs} base psi device
  SOURCES diago_david_float_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp 
          ../../module_basis/module_pw/test/test_tool.cpp
          ../../module_hamilt_general/operator.cpp
          ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
)
if(ENABLE_LCAO)
AddTest(
  TARGET HSolver_dav_real
  LIBS ${math_libs} base psi device
  SOURCES diago_david_real_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp 
          ../../module_basis/module_pw/test/test_tool.cpp
          ../../module_hamilt_general/operator.cpp
          ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
)
endif()

AddTest(
  TARGET HSolver_base
  LIBS ${math_libs} psi device base
  SOURCES test_hsolver.cpp 
)

AddTest(
  TARGET HSolver_pw
  LIBS ${math_libs} psi device base container
  SOURCES test_hsolver_pw.cpp ../hsolver_pw.cpp ../diago_bpcg.cpp ../diago_iter_assist.cpp
)

AddTest(
  TARGET HSolver_sdft
  LIBS ${math_libs} psi device base container
  SOURCES test_hsolver_sdft.cpp ../hsolver_pw_sdft.cpp ../hsolver_pw.cpp ../diago_bpcg.cpp ../diago_iter_assist.cpp
)

if(ENABLE_LCAO)
  if(USE_ELPA)
  AddTest(
    TARGET HSolver_LCAO
    LIBS ${math_libs} ELPA::ELPA base genelpa psi device
    SOURCES diago_lcao_test.cpp ../diago_elpa.cpp ../diago_blas.cpp 
  )
  else()
    AddTest(
      TARGET HSolver_LCAO
      LIBS ${math_libs} base psi device
      SOURCES diago_lcao_test.cpp ../diago_blas.cpp 
    )
  endif()
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