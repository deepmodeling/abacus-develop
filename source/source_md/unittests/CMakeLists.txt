abacus_disable_feature_definitions(__MPI __LCAO)
abacus_add_local_feature_definitions(__NORMAL)

addtest(
  TARGET MODULE_MD_LJ_potential
  LIBS parameter psi device
  SOURCES test_lj.cpp
          ../../source_cell/module_neighlist/test/test_domain_decomposition.cpp
          ${md_test_depend_files})

addtest(
  TARGET MODULE_MD_SNAP_potential
  LIBS parameter psi device
  SOURCES test_snap_potential.cpp
          ../module_snap/snap_potential.cpp
          ../module_snap/snap_calculator.cpp
          ../module_snap/snap_descriptor.cpp
          ../module_snap/snap_derivative.cpp
          ../module_snap/snap_model.cpp
          ../module_snap/snap_cg.cpp
          ${md_test_depend_files})

set(_md_dp_test_libs parameter psi device)
if(DEFINED DeePMD_DIR)
  if(DeePMDC_FOUND)
    list(APPEND _md_dp_test_libs DeePMD::deepmd_c)
  else()
    list(APPEND _md_dp_test_libs DeePMD::deepmd_cc)
  endif()
endif()

addtest(
  TARGET MODULE_MD_DP_potential
  LIBS ${_md_dp_test_libs}
  SOURCES test_dp.cpp ../module_dp/dp.cpp
          ../../source_cell/module_neighlist/neighbor_list.cpp
          ${md_test_depend_files})

set(_md_nep_test_libs parameter psi device)
if(DEFINED NEP_DIR AND NEP_FOUND)
  list(APPEND _md_nep_test_libs NEP::nep)
endif()

addtest(
  TARGET MODULE_MD_NEP_potential
  LIBS ${_md_nep_test_libs}
  SOURCES test_nep_potential.cpp ../module_nep/nep_potential.cpp
          ../../source_cell/module_neighlist/neighbor_list.cpp
          ${md_test_depend_files})

addtest(
  TARGET MODULE_MD_fire
  LIBS parameter psi device
  SOURCES test_fire.cpp ../md_base.cpp ../fire.cpp ${md_test_depend_files})

addtest(
  TARGET MODULE_MD_verlet
  LIBS parameter psi device
  SOURCES test_verlet.cpp ../test/test_md_base.cpp ../md_base.cpp ../verlet.cpp
          ${md_test_depend_files})

addtest(
  TARGET MODULE_MD_nhc
  LIBS parameter psi device
  SOURCES test_nhchain.cpp ../md_base.cpp ../nhchain.cpp
          ../../source_base/output.cpp ${md_test_depend_files})

addtest(
  TARGET MODULE_MD_msst
  LIBS parameter psi device
  SOURCES test_msst.cpp ../md_base.cpp ../msst.cpp
          ../../source_cell/update_cell.cpp ../../source_base/output.cpp
          ${md_test_depend_files})

addtest(
  TARGET MODULE_MD_lgv
  LIBS parameter psi device
  SOURCES test_langevin.cpp ../md_base.cpp ../langevin.cpp ${md_test_depend_files})

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_MD_dump_mpi
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter psi device
    SOURCES test_md_func.cpp ../../source_base/parallel_common.cpp
            ../../source_cell/magnetism.cpp ${md_test_depend_files})
  add_test(
    NAME MODULE_MD_dump_mpi_np4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_MD_dump_mpi> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  addtest(
    TARGET MODULE_MD_base_mpi
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter psi device
    SOURCES test_md_base.cpp ../md_base.cpp
            ../../source_base/parallel_common.cpp
            ../../source_cell/magnetism.cpp ${md_test_depend_files})
  add_test(
    NAME MODULE_MD_base_mpi_np4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_MD_base_mpi> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(
    MODULE_MD_dump_mpi MODULE_MD_dump_mpi_np4 MODULE_MD_base_mpi
    MODULE_MD_base_mpi_np4
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;I_MPI_PIN=0"
               TIMEOUT 60)
  foreach(rank_count IN ITEMS 1 4)
    math(EXPR processor_count "${rank_count} * 4")
    add_test(
      NAME MODULE_MD_base_mpi_np${rank_count}_omp4
      COMMAND
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank_count}
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:MODULE_MD_base_mpi>
        ${MPIEXEC_POSTFLAGS}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(
      MODULE_MD_base_mpi_np${rank_count}_omp4
      PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4;MKL_NUM_THREADS=1;I_MPI_PIN=0"
                 PROCESSORS "${processor_count}" TIMEOUT 60)
  endforeach()
endif()

# Separate MD -> RT contract; the MD unit target above keeps no RT sources.
# Register after MD support sources are available; the test lives with RT.
if(ENABLE_LCAO AND ENABLE_MPI)
  addtest(
    TARGET INTEGRATION_MD_RT_boundary
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter psi device container
    SOURCES ../../source_lcao/module_rt/test/test_md_rt_boundary.cpp
            ../../source_base/test/mpi_test_main.cpp
            ../md_base.cpp
            ../../source_lcao/module_rt/boundary_fix.cpp
            ../../source_lcao/module_rt/td_folding.cpp
            ../../source_hamilt/module_hcontainer/base_matrix.cpp
            ../../source_hamilt/module_hcontainer/atom_pair.cpp
            ../../source_hamilt/module_hcontainer/hcontainer.cpp
            ../../source_basis/module_ao/parallel_orbitals.cpp
            ../../source_base/block_cyclic_distribution.cpp
            ../../source_cell/magnetism.cpp
            ../../source_base/parallel_common.cpp
            ${md_test_depend_files})
  # One executable, separate CTest entries for the RT unit and MD/RT contract.
  abacus_add_rt_test_variants(
    INTEGRATION_MD_RT_boundary INTEGRATION_MD_RT_boundary "integration;md;rt"
    "*MdRtBoundaryTest.*")
  abacus_add_rt_test_variants(
    MODULE_LCAO_tddft_td_folding INTEGRATION_MD_RT_boundary "unit;rt"
    "TdFoldingTest.*")
endif()
