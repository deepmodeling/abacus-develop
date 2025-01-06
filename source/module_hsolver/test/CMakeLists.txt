remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-DUSE_PAW)
remove_definitions(-D__EXX)

if (ENABLE_MPI)
  AddTest(
    TARGET hsolver_parak2d_test
    LIBS parameter  ${math_libs} base device MPI::MPI_CXX
    SOURCES parallel_k2d_test.cpp ../parallel_k2d.cpp ../../module_cell/parallel_kpoints.cpp
  )
  AddTest(
    TARGET HSolver_bpcg
    LIBS parameter  ${math_libs} base psi device container
    SOURCES diago_bpcg_test.cpp ../diago_bpcg.cpp  ../diago_iter_assist.cpp  
            ../../module_basis/module_pw/test/test_tool.cpp
            ../../module_hamilt_general/operator.cpp
            ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
  )
  AddTest(
    TARGET HSolver_cg
    LIBS parameter  ${math_libs} base psi device container
    SOURCES diago_cg_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp  ../diag_const_nums.cpp
            ../../module_basis/module_pw/test/test_tool.cpp
            ../../module_hamilt_general/operator.cpp
            ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
  )
  AddTest(
    TARGET HSolver_cg_float
    LIBS parameter  ${math_libs} base psi device container
    SOURCES diago_cg_float_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp   ../diag_const_nums.cpp
            ../../module_basis/module_pw/test/test_tool.cpp
            ../../module_hamilt_general/operator.cpp
            ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
  )
  AddTest(
    TARGET HSolver_dav
    LIBS parameter  ${math_libs} base psi device
    SOURCES diago_david_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp  ../diag_const_nums.cpp
            ../../module_basis/module_pw/test/test_tool.cpp
            ../../module_hamilt_general/operator.cpp
            ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
  )
  AddTest(
    TARGET HSolver_dav_float
    LIBS parameter  ${math_libs} base psi device
    SOURCES diago_david_float_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp  ../diag_const_nums.cpp
            ../../module_basis/module_pw/test/test_tool.cpp
            ../../module_hamilt_general/operator.cpp
            ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
  )
  if(ENABLE_LCAO)
  AddTest(
    TARGET HSolver_cg_real
    LIBS parameter  ${math_libs} base psi device container
    SOURCES diago_cg_float_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp   ../diag_const_nums.cpp
            ../../module_basis/module_pw/test/test_tool.cpp
            ../../module_hamilt_general/operator.cpp
            ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
  )
  AddTest(
    TARGET HSolver_dav_real
    LIBS parameter  ${math_libs} base psi device
    SOURCES diago_david_real_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp  ../diag_const_nums.cpp
            ../../module_basis/module_pw/test/test_tool.cpp
            ../../module_hamilt_general/operator.cpp
            ../../module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.cpp
  )
  endif()

  AddTest(
    TARGET HSolver_base
    LIBS parameter  ${math_libs} psi device base
    SOURCES test_hsolver.cpp 
  )

  AddTest(
    TARGET HSolver_pw
    LIBS parameter  ${math_libs} psi device base container
    SOURCES test_hsolver_pw.cpp ../hsolver_pw.cpp ../hsolver_lcaopw.cpp ../diago_bpcg.cpp ../diago_dav_subspace.cpp ../diag_const_nums.cpp ../diago_iter_assist.cpp
  )

  AddTest(
    TARGET HSolver_sdft
    LIBS parameter  ${math_libs} psi device base container
    SOURCES test_hsolver_sdft.cpp ../hsolver_pw_sdft.cpp ../hsolver_pw.cpp ../diago_bpcg.cpp ../diago_dav_subspace.cpp ../diag_const_nums.cpp ../diago_iter_assist.cpp
  )

  if(ENABLE_LCAO)
    if(USE_ELPA)
    AddTest(
      TARGET HSolver_LCAO
      LIBS parameter  ${math_libs} ELPA::ELPA base genelpa psi device
      SOURCES diago_lcao_test.cpp ../diago_elpa.cpp ../diago_scalapack.cpp 
    )
    else()
      AddTest(
        TARGET HSolver_LCAO
        LIBS parameter  ${math_libs} base psi device
        SOURCES diago_lcao_test.cpp ../diago_scalapack.cpp 
      )
    endif()

    if (ENABLE_PEXSI)
    AddTest(
      TARGET HSolver_LCAO_PEXSI
      LIBS parameter  ${math_libs} ${PEXSI_LIBRARY} ${SuperLU_DIST_LIBRARY} ${ParMETIS_LIBRARY} ${METIS_LIBRARY} MPI::MPI_CXX base psi device pexsi
      SOURCES diago_pexsi_test.cpp ../diago_pexsi.cpp ../../module_basis/module_ao/parallel_orbitals.cpp
    )
    endif()
  endif()
  if (USE_CUDA)
  AddTest(
    TARGET HSolver_LCAO_cusolver
    LIBS parameter  ${math_libs}  base psi device 
    SOURCES diago_lcao_cusolver_test.cpp ../diago_cusolver.cpp ../diago_scalapack.cpp 
    ../kernels/math_kernel_op.cpp
    ../kernels/dngvd_op.cpp
    ../kernels/cuda/diag_cusolver.cu
  )
  endif()
else()
  if(ENABLE_LCAO)
  AddTest(
        TARGET HSolver_Lapack
        LIBS parameter  ${math_libs} base psi device
        SOURCES diago_lapack_test.cpp ../diago_lapack.cpp 
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

install(FILES GammaOnly-Si2-Solution.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES GammaOnly-Si64-Solution.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES KPoints-Si2-Solution.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES KPoints-Si64-Solution.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(FILES diago_cg_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES diago_david_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES diago_lcao_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(FILES PEXSI-H-GammaOnly-Si2.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES PEXSI-S-GammaOnly-Si2.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES PEXSI-DM-GammaOnly-Si2.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES diago_pexsi_parallel_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES parallel_k2d_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})



AddTest(
  TARGET hsolver_diago_hs_parallel
  LIBS parameter  ${math_libs} ELPA::ELPA base device MPI::MPI_CXX genelpa psi
  SOURCES test_diago_hs_para.cpp ../diag_hs_para.cpp ../diago_pxxxgvx.cpp ../diago_elpa.cpp ../diago_scalapack.cpp 
)

find_program(BASH bash)
if (ENABLE_MPI)
  add_test(NAME HSolver_cg_parallel
        COMMAND ${BASH} diago_cg_parallel_test.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
  add_test(NAME hsolver_parak2d_test_parallel
        COMMAND ${BASH} parallel_k2d_test.sh
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
    if(ENABLE_PEXSI)
      add_test(NAME HSolver_LCAO_PEXSI_parallel
            COMMAND ${BASH} diago_pexsi_parallel_test.sh
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
    endif()
  endif()
endif()