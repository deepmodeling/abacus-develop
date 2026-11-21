abacus_disable_feature_definitions(__MPI __LCAO)
abacus_add_local_feature_definitions(__NORMAL)

list(
  APPEND
  depend_files
  ../md_func.cpp
  ../../source_cell/basecell.cpp
  ../../source_cell/unitcell.cpp
  ../../source_cell/unitcell_index.cpp
  ../../source_cell/unitcell_stats.cpp
  ../../source_cell/unitcell_setup.cpp
  ../../source_cell/update_cell.cpp
  ../../source_cell/bcast_cell.cpp
  ../../source_cell/print_cell.cpp
  ../../source_cell/atom_spec.cpp
  ../../source_cell/read_stru.cpp
  ../../source_cell/read_atom_species.cpp
  ../../source_cell/atom_pseudo.cpp
  ../../source_cell/read_atoms.cpp
  ../../source_cell/read_atoms_helper.cpp
  ../../source_cell/cell_tools.cpp
  ../../source_cell/pseudo.cpp
  ../../source_cell/read_pp.cpp
  ../../source_cell/read_pp_complete.cpp
  ../../source_cell/read_pp_upf100.cpp
  ../../source_cell/read_pp_upf201.cpp
  ../../source_cell/read_pp_vwr.cpp
  ../../source_cell/read_pp_blps.cpp
  ../../source_base/matrix3.cpp
  ../../source_base/matrix.cpp
  ../../source_base/timer.cpp
  ../../source_base/module_external/blas_connector_base.cpp
  ../../source_base/module_external/blas_connector_matrix.cpp
  ../../source_base/module_external/blas_connector_vector.cpp
  ../../source_base/memory_recorder.cpp
  ../../source_base/global_variable.cpp
  ../../source_base/global_function.cpp
  ../../source_base/global_file.cpp
  ../../source_base/tool_title.cpp
  ../../source_base/tool_check.cpp
  ../../source_base/tool_quit.cpp
  ../../source_base/intarray.cpp
  ../../source_base/realarray.cpp
  ../../source_base/complexarray.cpp
  ../../source_base/complexmatrix.cpp
  ../../source_base/global_variable.cpp
  ../../source_base/libm/branred.cpp
  ../../source_base/libm/sincos.cpp
  ../../source_base/math_integral.cpp
  ../../source_cell/module_neighbor/sltk_atom_arrange.cpp
  ../../source_cell/module_neighbor/sltk_atom.cpp
  ../../source_cell/module_neighbor/sltk_grid.cpp
  ../../source_cell/module_neighbor/sltk_grid_driver.cpp
  ../../source_cell/module_neighlist/neighbor_search.cpp
  ../../source_cell/module_neighlist/bin_manager.cpp
  ../../source_cell/module_neighlist/page_allocator.cpp
  ../../source_cell/dom_core.cpp
  ../../source_cell/dom_ghost.cpp
  ../../source_cell/dom_migrate.cpp
  ../../source_cell/dom_force.cpp
  ../../source_cell/mdcell.cpp
  ../../source_cell/mdcell_adapter.cpp
  ../../source_cell/mdcell_reader.cpp
  ../../source_base/output.cpp
  ../../source_io/module_output/output_log.cpp
  ../../source_io/module_output/print_info.cpp
  ../../source_cell/cif_io.cpp
  ../../source_esolver/esolver_forcefields.cpp
  ../module_lj/lj_forcefield.cpp
  ../../source_base/parallel_reduce.cpp
  ../../source_base/parallel_global.cpp
  ../../source_base/parallel_cell.cpp
  ../../source_base/parallel_comm.cpp
  ../../source_cell/read_pp_ucell.cpp
  ../../source_cell/cal_wfc.cpp
  ../../source_cell/cal_nelec_nband.cpp
  ../../source_cell/read_orb.cpp
  ../../source_cell/sep.cpp
  ../../source_cell/sep_cell.cpp)

addtest(
  TARGET MODULE_MD_func
  LIBS parameter psi device
  SOURCES test_md_func.cpp ${depend_files})

addtest(
  TARGET MODULE_MD_SNAP_model
  SOURCES test_snap_model.cpp ../module_snap/snap_model.cpp)

addtest(
  TARGET MODULE_MD_SNAP_cg
  SOURCES test_snap_cg.cpp ../module_snap/snap_cg.cpp)

addtest(
  TARGET MODULE_MD_SNAP_descriptor
  LIBS base device
  SOURCES test_snap_descriptor.cpp ../module_snap/snap_descriptor.cpp
          ../module_snap/snap_derivative.cpp ../module_snap/snap_model.cpp
          ../module_snap/snap_cg.cpp)

addtest(
  TARGET MODULE_MD_SNAP_calculator
  LIBS base device
  SOURCES test_snap_calculator.cpp ../module_snap/snap_calculator.cpp
          ../module_snap/snap_descriptor.cpp ../module_snap/snap_derivative.cpp
          ../module_snap/snap_model.cpp ../module_snap/snap_cg.cpp)

addtest(
  TARGET MODULE_MD_lj_forcefield
  LIBS parameter psi device
  SOURCES test_lj_forcefield.cpp ${depend_files})

addtest(
  TARGET MODULE_MD_snap_forcefield
  LIBS parameter psi device
  SOURCES test_snap_forcefield.cpp
          ../module_snap/snap_forcefield.cpp
          ../module_snap/snap_calculator.cpp
          ../module_snap/snap_descriptor.cpp
          ../module_snap/snap_derivative.cpp
          ../module_snap/snap_model.cpp
          ../module_snap/snap_cg.cpp
          ${depend_files})

set(_md_dp_test_libs parameter psi device)
if(DEFINED DeePMD_DIR)
  if(DeePMDC_FOUND)
    list(APPEND _md_dp_test_libs DeePMD::deepmd_c)
  else()
    list(APPEND _md_dp_test_libs DeePMD::deepmd_cc)
  endif()
endif()

addtest(
  TARGET MODULE_MD_dp_forcefield
  LIBS ${_md_dp_test_libs}
  SOURCES test_dp_forcefield.cpp ../module_dp/dp_forcefield.cpp
          ../module_dp/dp_aparam.cpp
          ../../source_cell/module_neighlist/neighbor_list.cpp
          ${depend_files})

addtest(
  TARGET MODULE_MD_dp_aparam
  SOURCES test_dp_aparam.cpp ../module_dp/dp_aparam.cpp)

set(_md_nep_test_libs parameter psi device)
if(DEFINED NEP_DIR AND NEP_FOUND)
  list(APPEND _md_nep_test_libs NEP::nep)
endif()

addtest(
  TARGET MODULE_MD_nep_forcefield
  LIBS ${_md_nep_test_libs}
  SOURCES test_nep_forcefield.cpp ../module_nep/nep_forcefield.cpp
          ../../source_cell/module_neighlist/neighbor_list.cpp
          ${depend_files})

addtest(
  TARGET MODULE_MD_fire
  LIBS parameter psi device
  SOURCES test_fire.cpp ../md_base.cpp ../fire.cpp ${depend_files})

addtest(
  TARGET MODULE_MD_base
  LIBS parameter psi device
  SOURCES test_md_base.cpp ../md_base.cpp
          ../../source_cell/magnetism.cpp ${depend_files})

addtest(
  TARGET MODULE_MD_verlet
  LIBS parameter psi device
  SOURCES test_verlet.cpp ../md_base.cpp ../verlet.cpp
          ${depend_files})

addtest(
  TARGET MODULE_MD_nhc
  LIBS parameter psi device
  SOURCES test_nhchain.cpp ../md_base.cpp ../nhchain.cpp
          ../../source_base/output.cpp ${depend_files})

addtest(
  TARGET MODULE_MD_msst
  LIBS parameter psi device
  SOURCES test_msst.cpp ../md_base.cpp ../msst.cpp
          ../../source_cell/update_cell.cpp ../../source_base/output.cpp
          ${depend_files})

addtest(
  TARGET MODULE_MD_lgv
  LIBS parameter psi device
  SOURCES test_langevin.cpp ../md_base.cpp ../langevin.cpp ${depend_files})

# Header-only MD loop predicates (is_last_step / is_frame_write_step). Kept in
# its own translation unit so this target links without the heavy MD/cell stack.
addtest(
  TARGET MODULE_MD_final_step
  SOURCES test_run_md_final_step.cpp)

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_MD_lj_forcefield_mpi
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter psi device
    SOURCES test_lj_forcefield.cpp ../../source_base/test/mpi_test_main.cpp
            ../../source_base/parallel_common.cpp
            ${depend_files})
  set_tests_properties(
    MODULE_MD_lj_forcefield_mpi
    PROPERTIES ENVIRONMENT "GTEST_FILTER=LJPotential.SplitCommunicatorReducesWithinCell;OMP_NUM_THREADS=1"
               TIMEOUT 60)
  add_test(
    NAME MODULE_MD_lj_forcefield_mpi_np4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_MD_lj_forcefield_mpi>
            --gtest_filter=LJPotential.SplitCommunicatorReducesWithinCell
            ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(MODULE_MD_lj_forcefield_mpi_np4
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;I_MPI_PIN=0"
                                  TIMEOUT 60)
  addtest(
    TARGET MODULE_MD_dump_mpi
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter psi device
    SOURCES test_md_func.cpp ../../source_base/parallel_common.cpp
            ../../source_cell/magnetism.cpp ${depend_files})
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
            ../../source_cell/magnetism.cpp ${depend_files})
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
            ../../source_base/atom_basis_view.cpp
            ../../source_base/basis_flattened_indexing.cpp
            ../../source_base/parallel_blkcyc.cpp
            ../../source_base/element_basis_index.cpp
            ../../source_cell/magnetism.cpp
            ../../source_base/parallel_common.cpp
            ${depend_files})
  # One executable, separate CTest entries for the RT unit and MD/RT contract.
  abacus_add_rt_test_variants(
    INTEGRATION_MD_RT_boundary INTEGRATION_MD_RT_boundary "integration;md;rt"
    "*MdRtBoundaryTest.*")
  abacus_add_rt_test_variants(
    MODULE_LCAO_tddft_td_folding INTEGRATION_MD_RT_boundary "unit;rt"
    "TdFoldingTest.*")
endif()
