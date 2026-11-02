abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_HSOLVER_parak2d_test
    LIBS parameter base device MPI::MPI_CXX
    SOURCES test_parallel_k2d.cpp ../parallel_k2d.cpp
            ../../source_cell/parallel_kpoints.cpp)
  addtest(
    TARGET MODULE_HSOLVER_bpcg
    LIBS parameter base psi device container
    SOURCES test_diago_bpcg.cpp ../diago_bpcg.cpp ../para_lin_tf.cpp
            ../diago_iter_assist.cpp
            ../../source_basis/module_pw/test/test_tool.cpp)
  addtest(
    TARGET MODULE_HSOLVER_cg
    LIBS parameter base psi device container
    SOURCES test_diago_cg.cpp ../diago_cg.cpp ../diago_iter_assist.cpp
            ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp)
  addtest(
    TARGET MODULE_HSOLVER_cg_float
    LIBS parameter base psi device container
    SOURCES test_diago_cg_float.cpp ../diago_cg.cpp ../diago_iter_assist.cpp
            ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp)
  addtest(
    TARGET MODULE_HSOLVER_dav
    LIBS parameter base psi device
    SOURCES test_diago_david_double.cpp ../diago_david.cpp
            ../diago_iter_assist.cpp ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp)
  addtest(
    TARGET MODULE_HSOLVER_dav_float
    LIBS parameter base psi device
    SOURCES test_diago_david_float.cpp ../diago_david.cpp
            ../diago_iter_assist.cpp ../diag_const_nums.cpp
            ../../source_basis/module_pw/test/test_tool.cpp)
  if(ENABLE_LCAO)
    addtest(
      TARGET MODULE_HSOLVER_cg_real
      LIBS parameter base psi device container
      SOURCES test_diago_cg_float.cpp ../diago_cg.cpp ../diago_iter_assist.cpp
              ../diag_const_nums.cpp
              ../../source_basis/module_pw/test/test_tool.cpp)
    addtest(
      TARGET MODULE_HSOLVER_dav_real
      LIBS parameter base psi device
      SOURCES test_diago_david_real.cpp ../diago_david.cpp
              ../diago_iter_assist.cpp ../diag_const_nums.cpp
              ../../source_basis/module_pw/test/test_tool.cpp)
  endif()

  addtest(
    TARGET MODULE_HSOLVER_base
    LIBS parameter psi device base
    SOURCES test_hsolver.cpp ../hsolver.cpp)

  addtest(
    TARGET MODULE_HSOLVER_pw
    LIBS parameter psi device base container symmetry
    SOURCES test_hsolver_pw.cpp
            ../hsolver_pw.cpp
            ../hsolver_lcaopw.cpp
            ../diago_bpcg.cpp
            ../diago_dav_subspace.cpp
            ../diag_const_nums.cpp
            ../diago_iter_assist.cpp
            ../para_lin_tf.cpp
            ../../source_estate/elecstate_tools.cpp
            ../../source_estate/occupy.cpp
            ../../source_base/module_fft/fft_bundle.cpp
            ../../source_base/module_fft/fft_cpu.cpp
            ../../source_cell/klist.cpp
            ../../source_cell/klist_io.cpp
            ../../source_cell/parallel_kpoints.cpp
            ../../source_cell/reciprocal_grid.cpp)

  if(ENABLE_LCAO)
    if(TARGET ELPA::ELPA)
      addtest(
        TARGET MODULE_HSOLVER_LCAO
        LIBS parameter ELPA::ELPA base genelpa psi device
        SOURCES test_diago_lcao.cpp ../diago_elpa.cpp ../diago_elpa_native.cpp
                ../diago_scalapack.cpp ../diago_lapack.cpp)
    else()
      addtest(
        TARGET MODULE_HSOLVER_LCAO
        LIBS parameter base psi device
        SOURCES test_diago_lcao.cpp ../diago_scalapack.cpp ../diago_lapack.cpp)
    endif()

    if(TARGET PEXSI::PEXSI)
      addtest(
        TARGET MODULE_HSOLVER_LCAO_PEXSI
        LIBS parameter PEXSI::PEXSI base psi device pexsi
        SOURCES test_diago_pexsi.cpp ../diago_pexsi.cpp
                ../../source_basis/module_ao/parallel_orbitals.cpp)
    endif()
  endif()
  if(USE_CUDA)
    addtest(
      TARGET MODULE_HSOLVER_LCAO_cusolver
      LIBS parameter base psi device
      SOURCES test_diago_cusolver.cpp ../diago_cusolver.cpp
              ../diago_scalapack.cpp ../kernels/hegvd_op.cpp
              ../kernels/cuda/diag_cusolver.cu)
    target_compile_definitions(MODULE_HSOLVER_LCAO_cusolver PRIVATE __CUDA)
  endif()
endif()
abacus_copy_test_files(
  H-KPoints-Si2.dat
  H-GammaOnly-Si2.dat
  S-KPoints-Si2.dat
  S-GammaOnly-Si2.dat
  H-KPoints-Si64.dat
  H-GammaOnly-Si64.dat
  S-KPoints-Si64.dat
  S-GammaOnly-Si64.dat
  GammaOnly-Si2-Solution.dat
  GammaOnly-Si64-Solution.dat
  KPoints-Si2-Solution.dat
  KPoints-Si64-Solution.dat
  diago_cg_parallel_test.sh
  diago_david_parallel_test.sh
  diago_lcao_parallel_test.sh
  PEXSI-H-GammaOnly-Si2.dat
  PEXSI-S-GammaOnly-Si2.dat
  PEXSI-DM-GammaOnly-Si2.dat
  diago_pexsi_parallel_test.sh
  parallel_k2d_test.sh)

if(ENABLE_MPI)
  if(TARGET ELPA::ELPA)
    addtest(
      TARGET MODULE_HSOLVER_diago_hs_parallel
      LIBS parameter
           ELPA::ELPA
           base
           device
           MPI::MPI_CXX
           genelpa
           psi
      SOURCES test_diag_hs_para.cpp ../diag_hs_para.cpp ../diago_pxxxgvx.cpp
              ../diago_elpa.cpp ../diago_scalapack.cpp)
  else()
    addtest(
      TARGET MODULE_HSOLVER_diago_hs_parallel
      LIBS parameter base device MPI::MPI_CXX psi
      SOURCES test_diag_hs_para.cpp ../diag_hs_para.cpp ../diago_pxxxgvx.cpp
              ../diago_scalapack.cpp)
  endif()

  addtest(
    TARGET MODULE_HSOLVER_linear_trans
    LIBS parameter base device MPI::MPI_CXX
    SOURCES test_para_lin_tf.cpp ../para_lin_tf.cpp)

  add_test(
    NAME MODULE_HSOLVER_para_linear_trans
    COMMAND mpirun -np 4 "$<TARGET_FILE:MODULE_HSOLVER_linear_trans>"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

find_program(BASH bash)
if(ENABLE_MPI)
  add_test(
    NAME MODULE_HSOLVER_cg_parallel
    COMMAND ${BASH} diago_cg_parallel_test.sh "$<TARGET_FILE:MODULE_HSOLVER_cg>"
            "$<TARGET_FILE:MODULE_HSOLVER_cg_float>"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  add_test(
    NAME MODULE_HSOLVER_parak2d_test_parallel
    COMMAND ${BASH} parallel_k2d_test.sh
            "$<TARGET_FILE:MODULE_HSOLVER_parak2d_test>"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  add_test(
    NAME MODULE_HSOLVER_dav_parallel
    COMMAND
      ${BASH} diago_david_parallel_test.sh "$<TARGET_FILE:MODULE_HSOLVER_dav>"
      "$<TARGET_FILE:MODULE_HSOLVER_dav_float>"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  if(ENABLE_LCAO)
    add_test(
      NAME MODULE_HSOLVER_LCAO_parallel
      COMMAND
        ${BASH} diago_lcao_parallel_test.sh "$<TARGET_FILE:MODULE_HSOLVER_LCAO>"
        "$<TARGET_FILE:MODULE_HSOLVER_dav_real>"
        "$<TARGET_FILE:MODULE_HSOLVER_cg_real>"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    if(ENABLE_PEXSI)
      add_test(
        NAME MODULE_HSOLVER_LCAO_PEXSI_parallel
        COMMAND ${BASH} diago_pexsi_parallel_test.sh
                "$<TARGET_FILE:MODULE_HSOLVER_LCAO_PEXSI>"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()
  endif()
endif()

if(ENABLE_LCAO)
  set(DAVIDSON_OWNERSHIP_SOURCES
      test_diago_david.cpp ../diago_david.cpp ../diago_dav_subspace.cpp
      ../diago_cg.cpp ../diago_iter_assist.cpp ../diag_const_nums.cpp)
  set(DAVIDSON_OWNERSHIP_LIBS parameter base psi device container)
  if(ENABLE_MPI)
    list(APPEND DAVIDSON_OWNERSHIP_SOURCES ../diag_hs_para.cpp
         ../diago_pxxxgvx.cpp ../diago_scalapack.cpp)
    if(TARGET ELPA::ELPA)
      list(APPEND DAVIDSON_OWNERSHIP_SOURCES ../diago_elpa.cpp)
      list(APPEND DAVIDSON_OWNERSHIP_LIBS ELPA::ELPA genelpa)
    endif()
  endif()
  addtest(
    TARGET MODULE_HSOLVER_empty_blocks
    LIBS ${DAVIDSON_OWNERSHIP_LIBS}
    SOURCES ${DAVIDSON_OWNERSHIP_SOURCES})
  set_tests_properties(
    MODULE_HSOLVER_empty_blocks
    PROPERTIES TIMEOUT 120 ENVIRONMENT
               "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  if(ENABLE_MPI)
    add_test(NAME MODULE_HSOLVER_empty_blocks_mpi
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                     $<TARGET_FILE:MODULE_HSOLVER_empty_blocks>)
    set_tests_properties(
      MODULE_HSOLVER_empty_blocks_mpi
      PROPERTIES PROCESSORS 4 TIMEOUT 120 ENVIRONMENT
                 "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  endif()
endif()
