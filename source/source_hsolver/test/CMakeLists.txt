abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

if (ENABLE_MPI)
  AddTest(
    TARGET MODULE_HSOLVER_parak2d_test
    LIBS parameter  base device MPI::MPI_CXX
    SOURCES parallel_k2d_test.cpp ../parallel_k2d.cpp ../../source_cell/parallel_kpoints.cpp
  )
  AddTest(
    TARGET MODULE_HSOLVER_bpcg
    LIBS parameter  base psi device container
    SOURCES diago_bpcg_test.cpp ../diago_bpcg.cpp ../para_linear_transform.cpp  ../diago_iter_assist.cpp  
            ../../source_basis/module_pw/test/test_tool.cpp
            ../../source_hamilt/operator.cpp
            ../../source_pw/module_pwdft/op_pw.cpp
  )
  AddTest(
    TARGET MODULE_HSOLVER_cg
    LIBS parameter  base psi device container
    SOURCES diago_cg_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp  ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp
            ../../source_hamilt/operator.cpp
            ../../source_pw/module_pwdft/op_pw.cpp
  )
  AddTest(
    TARGET MODULE_HSOLVER_cg_float
    LIBS parameter  base psi device container
    SOURCES diago_cg_float_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp   ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp
            ../../source_hamilt/operator.cpp
            ../../source_pw/module_pwdft/op_pw.cpp
  )
  AddTest(
    TARGET MODULE_HSOLVER_dav
    LIBS parameter  base psi device
    SOURCES diago_david_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp  ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp
            ../../source_hamilt/operator.cpp
            ../../source_pw/module_pwdft/op_pw.cpp
  )
  AddTest(
    TARGET MODULE_HSOLVER_dav_float
    LIBS parameter  base psi device
    SOURCES diago_david_float_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp  ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp
            ../../source_hamilt/operator.cpp
            ../../source_pw/module_pwdft/op_pw.cpp
  )
  if(ENABLE_LCAO)
  AddTest(
    TARGET MODULE_HSOLVER_cg_real
    LIBS parameter  base psi device container
    SOURCES diago_cg_float_test.cpp ../diago_cg.cpp  ../diago_iter_assist.cpp   ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp
            ../../source_hamilt/operator.cpp
            ../../source_pw/module_pwdft/op_pw.cpp
  )
  AddTest(
    TARGET MODULE_HSOLVER_dav_real
    LIBS parameter  base psi device
    SOURCES diago_david_real_test.cpp ../diago_david.cpp  ../diago_iter_assist.cpp  ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp
            ../../source_hamilt/operator.cpp
            ../../source_pw/module_pwdft/op_pw.cpp
  )
  endif()

  AddTest(
    TARGET MODULE_HSOLVER_base
    LIBS parameter  psi device base
    SOURCES test_hsolver.cpp 
  )

  AddTest(
    TARGET MODULE_HSOLVER_pw
    LIBS parameter  psi device base container
    SOURCES test_hsolver_pw.cpp ../hsolver_pw.cpp ../hsolver_lcaopw.cpp ../diago_bpcg.cpp ../diago_dav_subspace.cpp ../diag_const_nums.cpp ../diago_iter_assist.cpp ../para_linear_transform.cpp
    ../../source_estate/elecstate_tools.cpp ../../source_estate/occupy.cpp ../../source_base/module_fft/fft_bundle.cpp ../../source_base/module_fft/fft_cpu.cpp
  )

  AddTest(
    TARGET MODULE_HSOLVER_sdft
    LIBS parameter  psi device base container
    SOURCES test_hsolver_sdft.cpp ../hsolver_pw_sdft.cpp ../hsolver_pw.cpp ../diago_bpcg.cpp ../diago_dav_subspace.cpp ../diag_const_nums.cpp ../diago_iter_assist.cpp ../para_linear_transform.cpp
                ../../source_estate/elecstate_tools.cpp ../../source_estate/occupy.cpp ../../source_base/module_fft/fft_bundle.cpp ../../source_base/module_fft/fft_cpu.cpp
    )

  if(ENABLE_LCAO)
    if(TARGET ELPA::ELPA)
    AddTest(
      TARGET MODULE_HSOLVER_LCAO
      LIBS parameter  ELPA::ELPA base genelpa psi device
      SOURCES diago_lcao_test.cpp ../diago_elpa.cpp ../diago_scalapack.cpp ../diago_lapack.cpp
    )
    else()
      AddTest(
        TARGET MODULE_HSOLVER_LCAO
        LIBS parameter  base psi device
        SOURCES diago_lcao_test.cpp ../diago_scalapack.cpp ../diago_lapack.cpp
      )
    endif()

    if (ENABLE_PEXSI)
    AddTest(
      TARGET MODULE_HSOLVER_LCAO_PEXSI
      LIBS parameter  ${PEXSI_LIBRARY} ${SuperLU_DIST_LIBRARY} ${ParMETIS_LIBRARY} ${METIS_LIBRARY} MPI::MPI_CXX base psi device pexsi
      SOURCES diago_pexsi_test.cpp ../diago_pexsi.cpp ../../source_basis/module_ao/parallel_orbitals.cpp
    )
    endif()
  endif()
  if (USE_CUDA)
  AddTest(
    TARGET MODULE_HSOLVER_LCAO_cusolver
    LIBS parameter   base psi device
    SOURCES diago_lcao_cusolver_test.cpp ../diago_cusolver.cpp ../diago_scalapack.cpp
    ../kernels/hegvd_op.cpp
    ../kernels/cuda/diag_cusolver.cu
  )
  target_compile_definitions(MODULE_HSOLVER_LCAO_cusolver PRIVATE __CUDA)
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


if(ENABLE_MPI)
  if(TARGET ELPA::ELPA)
    AddTest(
      TARGET MODULE_HSOLVER_diago_hs_parallel
      LIBS parameter ELPA::ELPA base device MPI::MPI_CXX genelpa psi
      SOURCES
        test_diago_hs_para.cpp
        ../diag_hs_para.cpp
        ../diago_pxxxgvx.cpp
        ../diago_elpa.cpp
        ../diago_scalapack.cpp
    )
  else()
    AddTest(
      TARGET MODULE_HSOLVER_diago_hs_parallel
      LIBS parameter base device MPI::MPI_CXX psi
      SOURCES
        test_diago_hs_para.cpp
        ../diag_hs_para.cpp
        ../diago_pxxxgvx.cpp
        ../diago_scalapack.cpp
    )
  endif()

  AddTest(
    TARGET MODULE_HSOLVER_linear_trans
    LIBS parameter base device MPI::MPI_CXX
    SOURCES test_para_linear_trans.cpp ../para_linear_transform.cpp
  )

  add_test(
    NAME MODULE_HSOLVER_para_linear_trans
    COMMAND mpirun -np 4 ./MODULE_HSOLVER_linear_trans
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

find_program(BASH bash)
if (ENABLE_MPI)
  add_test(NAME MODULE_HSOLVER_cg_parallel
        COMMAND ${BASH} diago_cg_parallel_test.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
  add_test(NAME MODULE_HSOLVER_parak2d_test_parallel
        COMMAND ${BASH} parallel_k2d_test.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_test(NAME MODULE_HSOLVER_dav_parallel
        COMMAND ${BASH} diago_david_parallel_test.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      ) 
  if(ENABLE_LCAO)
    add_test(NAME MODULE_HSOLVER_LCAO_parallel
          COMMAND ${BASH} diago_lcao_parallel_test.sh
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
    if(ENABLE_PEXSI)
      add_test(NAME MODULE_HSOLVER_LCAO_PEXSI_parallel
            COMMAND ${BASH} diago_pexsi_parallel_test.sh
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
    endif()
  endif()
endif()
