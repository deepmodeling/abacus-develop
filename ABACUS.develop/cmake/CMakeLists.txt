########################################
# CMake build system
# This file is part of ABACUS
cmake_minimum_required(VERSION 3.18)
########################################

project(ABACUS
    VERSION 2.2.0
    DESCRIPTION "ABACUS is an electronic structure package based on DFT."
    HOMEPAGE_URL "https://github.com/deepmodeling/abacus-develop"
    LANGUAGES CXX
)
set(ABACUS_BIN_NAME ${PROJECT_NAME}-${PROJECT_VERSION})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

get_filename_component(ABACUS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)
set(ABACUS_SOURCE_DIR ${ABACUS_DIR}/source)

set(CEREAL_DIR /usr/local)
include_directories(${CEREAL_DIR}/include)

set(ELPA_DIR /usr/local)
include_directories(${ELPA_DIR}/include)

set(FFTW_DIR /usr/local)
include_directories(${FFTW_DIR}/include)

set(SCALAPACK_DIR /usr/local)
include_directories(${SCALAPACK_DIR}/include)

# OpenBlas is required by ScaLAPACK
find_package(BLAS REQUIRED)

# We only use the header-only part of Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
add_compile_options(${MPI_CXX_COMPILE_FLAGS})

# Look for pthread or equivalent
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

set(EXTERNAL_LIB_NAMES
    -lgfortran
    -lm
    -lelpa
    -lfftw3
    ${SCALAPACK_DIR}/lib/libscalapack.a
)

set(EXTERNAL_LIB_DIRS
    ${LAPACK_DIR}/lib
    ${FFTW_DIR}/lib
    ${ELPA_DIR}/lib
)

add_definitions(
    -D__EXX
    -D__FFTW3
    -D__FP
    -D__MPI
    -D__OPENMP
    -D__SELINV
    -DMETIS
    -DEXX_DM=3
    -DEXX_H_COMM=2
    -DTEST_EXX_LCAO=0
    -DTEST_EXX_RADIAL=1
    -DUSE_CEREAL_SERIALIZATION
)

set(CODE_CONF
    driver.cpp
    input_conv.cpp
    input.cpp
    run_lcao.cpp
    run_pw.cpp
    src_io/write_input.cpp
)

SET(CODE_PW
    src_io/berryphase.cpp
    src_io/chi0_hilbert.cpp
    src_io/chi0_standard.cpp
    src_io/epsilon0_pwscf.cpp
    src_io/epsilon0_vasp.cpp
    src_io/read_pp.cpp
    src_io/read_pp_upf100.cpp
    src_io/read_pp_upf201.cpp
    src_io/read_pp_vwr.cpp
    src_io/to_wannier90.cpp
    src_io/unk_overlap_pw.cpp
    src_ions/MD_basic.cpp
    src_ions/MD_fire.cpp
    src_ions/MD_func.cpp
    src_ions/MD_thermo.cpp
    src_pw/VL_in_pw.cpp
    src_pw/VNL_in_pw.cpp
    src_pw/atom_pseudo.cpp
    src_pw/electrons.cpp
    src_pw/exx_lip.cpp
    src_pw/forces.cpp
    src_pw/pseudo_nc.cpp
    src_pw/soc.cpp
    src_pw/sto_che.cpp
    src_pw/sto_elec.cpp
    src_pw/sto_hchi.cpp
    src_pw/sto_iter.cpp
    src_pw/sto_wf.cpp
    src_pw/stress_func_cc.cpp
    src_pw/stress_func_ewa.cpp
    src_pw/stress_func_gga.cpp
    src_pw/stress_func_har.cpp
    src_pw/stress_func_kin.cpp
    src_pw/stress_func_loc.cpp
    src_pw/stress_func_nl.cpp
    src_pw/stress_func_print.cpp
    src_pw/stress_pw.cpp
    src_pw/threshold_elec.cpp
    src_pw/unitcell_pseudo.cpp
    src_pw/vdwd2.cpp
    src_pw/vdwd2_parameters.cpp
    src_pw/vdwd3.cpp
    src_pw/vdwd3_parameters.cpp
    src_pw/xc_1.cpp
    src_pw/xc_2.cpp
    src_pw/xc_3.cpp
    src_pw/xc_functional.cpp
    src_pw/xc_type.cpp
)

set(CODE_TOOLS
    src_global/complexarray.cpp
    src_global/complexmatrix.cpp
    src_global/element_basis_index.cpp
    src_global/export.cpp
    src_global/global_file.cpp
    src_global/global_function.cpp
    src_global/global_variable.cpp
    src_global/intarray.cpp
    src_global/integral.cpp
    src_global/math_integral.cpp
    src_global/math_polyint.cpp
    src_global/math_sphbes.cpp
    src_global/math_ylmreal.cpp
    src_global/mathzone_add1.cpp
    src_global/mathzone.cpp
    src_global/matrix.cpp
    src_global/matrix3.cpp
    src_global/mymath1.cpp
    src_global/mymath3.cpp
    src_global/memory.cpp
    src_global/realarray.cpp
    src_global/restart.cpp
    src_global/timer.cpp
    src_global/poission.cpp
    src_global/polint.cpp
    src_global/sph_bessel_recursive-d1.cpp
    src_global/sph_bessel_recursive-d2.cpp
    src_global/sph_bessel.cpp
    src_global/tool_check.cpp
    src_global/tool_quit.cpp
    src_global/tool_title.cpp
    src_io/print_info.cpp
)

set(CODE_LCAO
    module_ORB/ORB_atomic.cpp
    module_ORB/ORB_atomic_lm.cpp
    module_ORB/ORB_control.cpp
    module_ORB/ORB_gaunt_table.cpp
    module_ORB/ORB_gen_tables.cpp
    module_ORB/ORB_nonlocal.cpp
    module_ORB/ORB_nonlocal_lm.cpp
    module_ORB/ORB_read.cpp
    module_ORB/ORB_table_alpha.cpp
    module_ORB/ORB_table_beta.cpp
    module_ORB/ORB_table_phi.cpp
    src_global/sltk_adjacent_set.cpp
    src_global/sltk_atom.cpp
    src_global/sltk_atom_arrange.cpp
    src_global/sltk_atom_input.cpp
    src_global/sltk_grid.cpp
    src_global/sltk_grid_driver.cpp
    src_global/ylm.cpp
    src_io/bessel_basis.cpp
    src_io/cal_r_overlap_R.cpp
    src_io/istate_charge.cpp
    src_io/istate_envelope.cpp
    src_io/numerical_basis.cpp
    src_io/numerical_descriptor.cpp
    src_io/unk_overlap_lcao.cpp
    src_lcao/DM_gamma.cpp
    src_lcao/DM_k.cpp
    src_lcao/ELEC_cbands_gamma.cpp
    src_lcao/ELEC_cbands_k.cpp
    src_lcao/ELEC_evolve.cpp
    src_lcao/ELEC_nscf.cpp
    src_lcao/ELEC_scf.cpp
    src_lcao/FORCE_STRESS.cpp
    src_lcao/FORCE_gamma.cpp
    src_lcao/FORCE_gamma_edm.cpp
    src_lcao/FORCE_gamma_tvnl.cpp
    src_lcao/FORCE_gamma_vl.cpp
    src_lcao/FORCE_k.cpp
    src_lcao/LCAO_descriptor.cpp
    src_lcao/LCAO_diago.cpp
    src_lcao/LCAO_evolve.cpp
    src_lcao/LCAO_gen_fixedH.cpp
    src_lcao/LCAO_hamilt.cpp
    src_lcao/LCAO_matrix.cpp
    src_lcao/LCAO_nnr.cpp
    src_lcao/LOOP_cell.cpp
    src_lcao/LOOP_elec.cpp
    src_lcao/LOOP_ions.cpp
    src_lcao/build_st_pw.cpp
    src_lcao/center2_orb-orb11.cpp
    src_lcao/center2_orb-orb21.cpp
    src_lcao/center2_orb-orb22.cpp
    src_lcao/gint_gamma.cpp
    src_lcao/gint_gamma_common.cpp
    src_lcao/gint_gamma_env.cpp
    src_lcao/gint_gamma_fvl.cpp
    src_lcao/gint_gamma_mull.cpp
    src_lcao/gint_gamma_rho.cpp
    src_lcao/gint_gamma_vl.cpp
    src_lcao/gint_k.cpp
    src_lcao/gint_k_fvl.cpp
    src_lcao/gint_k_init.cpp
    src_lcao/gint_k_rho.cpp
    src_lcao/gint_k_vl.cpp
    src_lcao/gint_tools.cpp
    src_lcao/global_fp.cpp
    src_lcao/grid_base.cpp
    src_lcao/grid_base_beta.cpp
    src_lcao/grid_bigcell.cpp
    src_lcao/grid_meshball.cpp
    src_lcao/grid_meshcell.cpp
    src_lcao/grid_meshk.cpp
    src_lcao/grid_technique.cpp
    src_lcao/local_orbital_charge.cpp
    src_lcao/local_orbital_wfc.cpp
    src_lcao/record_adj.cpp
    src_lcao/run_md.cpp
    src_lcao/wfc_dm_2d.cpp
    src_ri/abfs-vector3_order.cpp
    src_ri/abfs.cpp
    src_ri/conv_coulomb_pot.cpp
    src_ri/conv_coulomb_pot_k.cpp
    src_ri/exx_abfs-abfs_index.cpp
    src_ri/exx_abfs-construct_orbs.cpp
    src_ri/exx_abfs-dm.cpp
    src_ri/exx_abfs-inverse_matrix_double.cpp
    src_ri/exx_abfs-io.cpp
    src_ri/exx_abfs-jle.cpp
    src_ri/exx_abfs-matrix_lcaoslcaos_lcaoslcaos.cpp
    src_ri/exx_abfs-matrix_orbs11.cpp
    src_ri/exx_abfs-matrix_orbs21.cpp
    src_ri/exx_abfs-matrix_orbs22.cpp
    src_ri/exx_abfs-parallel-communicate-dm3-allreduce.cpp
    src_ri/exx_abfs-parallel-communicate-dm3.cpp
    src_ri/exx_abfs-parallel-communicate-function.cpp
    src_ri/exx_abfs-parallel-communicate-hexx-allreduce2.cpp
    src_ri/exx_abfs-parallel-communicate-hexx.cpp
    src_ri/exx_abfs-parallel-distribute-htime.cpp
    src_ri/exx_abfs-parallel-distribute-kmeans.cpp
    src_ri/exx_abfs-parallel-distribute-order.cpp
    src_ri/exx_abfs-pca.cpp
    src_ri/exx_abfs-screen-cauchy.cpp
    src_ri/exx_abfs-screen-schwarz.cpp
    src_ri/exx_abfs-util.cpp
    src_ri/exx_abfs.cpp
    src_ri/exx_lcao.cpp
    src_ri/exx_opt_orb-print.cpp
    src_ri/exx_opt_orb.cpp
    src_parallel/parallel_orbitals.cpp
    src_parallel/subgrid_oper.cpp
    src_pdiag/GenELPA.cpp
    src_pdiag/pdiag_basic.cpp
    src_pdiag/pdiag_common.cpp
    src_pdiag/pdiag_double.cpp
    src_pw/wavefunc_in_pw.cpp
)

set(CODE_PARALLEL
    src_parallel/parallel_global.cpp
    src_parallel/parallel_kpoints.cpp
    src_parallel/parallel_common.cpp
    src_parallel/parallel_reduce.cpp
    src_parallel/parallel_pw.cpp
    src_parallel/ft.cpp
    src_parallel/parallel_grid.cpp
)

set(CODE_FIRST_PRINCIPLES
    ${CODE_CONF}
    ${CODE_PW}
    ${CODE_LCAO}
    input_update.cpp
    src_io/cal_test.cpp
    src_io/optical.cpp
    src_io/read_atoms.cpp
    src_io/read_cell_pseudopots.cpp
    src_io/read_dm.cpp
    src_io/read_rho.cpp
    src_io/winput.cpp
    src_io/write_HS.cpp
    src_io/write_HS_R.cpp
    src_io/write_dm.cpp
    src_io/write_pot.cpp
    src_io/write_rho.cpp
    src_io/write_rho_dipole.cpp
    src_ions/bfgs_basic.cpp
    src_ions/ions_move_basic.cpp
    src_ions/ions_move_bfgs.cpp
    src_ions/ions_move_cg.cpp
    src_ions/ions_move_methods.cpp
    src_ions/ions_move_sd.cpp
    src_ions/ions.cpp
    src_ions/lattice_change_basic.cpp
    src_ions/lattice_change_cg.cpp
    src_ions/lattice_change_methods.cpp
    src_ions/variable_cell.cpp
    src_lcao/dftu.cpp
    src_lcao/dftu_relax.cpp
    src_lcao/dftu_yukawa.cpp
    src_pw/H_Ewald_pw.cpp
    src_pw/H_Hartree_pw.cpp
    src_pw/H_TDDFT_pw.cpp
    src_pw/H_XC_pw.cpp
    src_pw/charge.cpp
    src_pw/charge_broyden.cpp
    src_pw/charge_extra.cpp
    src_pw/charge_mixing.cpp
    src_pw/charge_pulay.cpp
    src_pw/efield.cpp
    src_pw/magnetism.cpp
    src_pw/occupy.cpp
    src_pw/potential.cpp
    src_pw/potential_libxc.cpp
    src_pw/xc_gga_pw.cpp
)

set(CODE_COMMON
    src_global/inverse_matrix.cpp
    src_io/dos.cpp
    src_io/energy_dos.cpp
    src_io/eximport.cpp
    src_io/mulliken_charge.cpp
    src_io/output.cpp
    src_io/rwstream.cpp
    src_io/wf_io.cpp
    src_io/wf_local.cpp
    src_pw/atom_spec.cpp
    src_pw/diago_cg.cpp
    src_pw/diago_david.cpp
    src_pw/energy.cpp
    src_pw/global.cpp
    src_pw/hamilt.cpp
    src_pw/hamilt_pw.cpp
    src_pw/klist.cpp
    src_pw/pw_basis.cpp
    src_pw/pw_complement.cpp
    src_pw/symm_other.cpp
    src_pw/symmetry.cpp
    src_pw/symmetry_basic.cpp
    src_pw/symmetry_rho.cpp
    src_pw/unitcell.cpp
    src_pw/use_fft.cpp
    src_pw/wavefunc.cpp
    src_pw/wf_atomic.cpp
    src_pw/wf_igk.cpp
)

set(CODE_PDIAG
    src_pdiag/pdgseps.cpp
    src_pdiag/pdst2g.cpp
    src_pdiag/pdstebz.cpp
    src_pdiag/pdsteiz.cpp
    src_pdiag/pdsyg2st.cpp
    src_pdiag/pdsytrd.cpp
    src_pdiag/pdt2s.cpp
    src_pdiag/pdtrsm.cpp
    src_pdiag/pzgseps.cpp
    src_pdiag/pzheg2st.cpp
    src_pdiag/pzhetrd.cpp
    src_pdiag/pzhtrsm.cpp
    src_pdiag/pzst2g.cpp
    src_pdiag/pzsteiz.cpp
    src_pdiag/pzt2s.cpp
)

set(CODE_PDIAG_MR
    src_pdiag/MRRR/dcopy.cpp
    src_pdiag/MRRR/dlae2.cpp
    src_pdiag/MRRR/dlaebz.cpp
    src_pdiag/MRRR/dlaev2.cpp
    src_pdiag/MRRR/dlaneg.cpp
    src_pdiag/MRRR/dlanst.cpp
    src_pdiag/MRRR/dlar1v.cpp
    src_pdiag/MRRR/dlarnv.cpp
    src_pdiag/MRRR/dlarra.cpp
    src_pdiag/MRRR/dlarrb.cpp
    src_pdiag/MRRR/dlarrc.cpp
    src_pdiag/MRRR/dlarrd.cpp
    src_pdiag/MRRR/dlarre.cpp
    src_pdiag/MRRR/dlarrf.cpp
    src_pdiag/MRRR/dlarrj.cpp
    src_pdiag/MRRR/dlarrk.cpp
    src_pdiag/MRRR/dlarrr.cpp
    src_pdiag/MRRR/dlarrv.cpp
    src_pdiag/MRRR/dlaruv.cpp
    src_pdiag/MRRR/dlas2.cpp
    src_pdiag/MRRR/dlascl.cpp
    src_pdiag/MRRR/dlaset.cpp
    src_pdiag/MRRR/dlasq2.cpp
    src_pdiag/MRRR/dlasq3.cpp
    src_pdiag/MRRR/dlasq4.cpp
    src_pdiag/MRRR/dlasq5.cpp
    src_pdiag/MRRR/dlasq6.cpp
    src_pdiag/MRRR/dlasrt.cpp
    src_pdiag/MRRR/dlassq.cpp
    src_pdiag/MRRR/dscal.cpp
    src_pdiag/MRRR/dstemr_mpi.cpp
    src_pdiag/MRRR/dswap.cpp
    src_pdiag/MRRR/i_nint.cpp
    src_pdiag/MRRR/ieeeck.cpp
    src_pdiag/MRRR/ilaenv.cpp
    src_pdiag/MRRR/iparmq.cpp
    src_pdiag/MRRR/lsame.cpp
    src_pdiag/MRRR/psort_w.cpp
    src_pdiag/MRRR/s_cmp.cpp
    src_pdiag/MRRR/s_copy.cpp
    src_pdiag/MRRR/xerbla.cpp
)

list(TRANSFORM CODE_COMMON PREPEND ${ABACUS_SOURCE_DIR}/)
list(TRANSFORM CODE_FIRST_PRINCIPLES PREPEND ${ABACUS_SOURCE_DIR}/)
list(TRANSFORM CODE_PARALLEL PREPEND ${ABACUS_SOURCE_DIR}/)
list(TRANSFORM CODE_TOOLS PREPEND ${ABACUS_SOURCE_DIR}/)
list(TRANSFORM CODE_PDIAG PREPEND ${ABACUS_SOURCE_DIR}/)
list(TRANSFORM CODE_PDIAG_MR PREPEND ${ABACUS_SOURCE_DIR}/)

add_library(abacusCommon OBJECT ${CODE_COMMON})
add_library(abacusFirstPrinciple OBJECT ${CODE_FIRST_PRINCIPLES})
add_library(abacusParallel OBJECT ${CODE_PARALLEL})
add_library(abacusTools OBJECT ${CODE_TOOLS})
add_library(abacusCodePdiag OBJECT ${CODE_PDIAG})
add_library(abacusCodePdiagMr OBJECT ${CODE_PDIAG_MR})

target_include_directories(abacusCommon PUBLIC ${ABACUS_SOURCE_DIR})
target_include_directories(abacusFirstPrinciple PUBLIC ${ABACUS_SOURCE_DIR})
target_include_directories(abacusParallel PUBLIC ${ABACUS_SOURCE_DIR})
target_include_directories(abacusTools PUBLIC ${ABACUS_SOURCE_DIR})
target_include_directories(abacusCodePdiag PUBLIC ${ABACUS_SOURCE_DIR})
target_include_directories(abacusCodePdiagMr PUBLIC ${ABACUS_SOURCE_DIR})

add_executable(${ABACUS_BIN_NAME}
    ${ABACUS_SOURCE_DIR}/main.cpp
    $<TARGET_OBJECTS:abacusCommon>
    $<TARGET_OBJECTS:abacusFirstPrinciple>
    $<TARGET_OBJECTS:abacusParallel>
    $<TARGET_OBJECTS:abacusTools>
    $<TARGET_OBJECTS:abacusCodePdiag>
    $<TARGET_OBJECTS:abacusCodePdiagMr>
)
target_link_libraries(${ABACUS_BIN_NAME}
    ${EXTERNAL_LIB_NAMES}
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    BLAS::BLAS
    Threads::Threads
)
target_link_directories(${ABACUS_BIN_NAME}
    PUBLIC ${EXTERNAL_LIB_DIRS}
)

