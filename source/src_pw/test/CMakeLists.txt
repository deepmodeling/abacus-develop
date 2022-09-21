remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

AddTest(
  TARGET abacus_pw.exe #Test PW ABACUS without LCAO #Please also update Makefile in src_pw
  LIBS ${math_libs} planewave base surchem xc_ operator_ks_pw neighbor symmetry
  SOURCES ../../driver.cpp ../../input.cpp ../../input_conv.cpp ../../driver_run.cpp ../../main.cpp ../../input_update.cpp
  # src_pw
  ../charge.cpp  ../charge_mixing.cpp  ../charge_pulay.cpp  ../charge_broyden.cpp  ../charge_extra.cpp  ../potential.cpp
  ../H_Ewald_pw.cpp  ../H_Hartree_pw.cpp     ../../src_io/write_wfc_realspace.cpp  ../magnetism.cpp    ../run_md_pw.cpp  ../occupy.cpp
    ../xc_3.cpp  ../vdwd2.cpp  ../vdwd2_parameters.cpp  ../vdwd3.cpp  ../vdwd3_parameters.cpp  ../VL_in_pw.cpp  ../VNL_in_pw.cpp  
    ../threshold_elec.cpp  ../electrons.cpp  ../forces.cpp  ../global.cpp  ../stress_pw.cpp  ../stress_func_cc.cpp  ../stress_func_gga.cpp
  ../stress_func_mgga.cpp  ../stress_func_kin.cpp  ../stress_func_nl.cpp  ../stress_func_ewa.cpp  ../stress_func_har.cpp  
  ../stress_func_loc.cpp  ../stress_func_print.cpp ../soc.cpp ../energy.cpp
    ../sto_iter.cpp  ../sto_wf.cpp  ../sto_func.cpp  ../sto_hchi.cpp
  ../sto_forces.cpp  ../sto_stress_pw.cpp ../structure_factor.cpp ../hamilt.cpp ../hamilt_pw_old.cpp ../pw_complement.cpp
  ../diago_david_old.cpp ../diago_cg_old.cpp ../wavefunc.cpp  ../wf_atomic.cpp ../wf_igk.cpp ../use_fft.cpp ../klist.cpp 
  ../symmetry_rho.cpp 
  # src_parallel
  ../../src_parallel/parallel_global.cpp ../../src_parallel/parallel_kpoints.cpp ../../src_parallel/parallel_common.cpp ../../src_parallel/parallel_reduce.cpp
  ../../src_parallel/parallel_grid.cpp
  # module_esolver
  ../../module_esolver/esolver.cpp  ../../module_esolver/esolver_ks.cpp ../../module_esolver/esolver_fp.cpp ../../module_esolver/esolver_ks_pw.cpp
  ../../module_esolver/esolver_ks_pw_tool.cpp ../../module_esolver/esolver_sdft_pw.cpp ../../module_esolver/esolver_sdft_pw_tool.cpp
  ../../module_esolver/esolver_lj.cpp  ../../module_esolver/esolver_dp.cpp
  # module_hsolver module_elecstate module_psi module_hamilt
  ../../module_hsolver/diago_cg.cpp ../../module_hsolver/diago_david.cpp ../../module_hsolver/diago_iter_assist.cpp
  ../../module_hsolver/hsolver_pw.cpp ../../module_hsolver/hsolver_pw_sdft.cpp
  ../../module_elecstate/elecstate.cpp ../../module_elecstate/elecstate_pw.cpp ../../module_elecstate/elecstate_pw_sdft.cpp
  ../../module_psi/psi.cpp
  ../../module_hamilt/hamilt_pw.cpp
  # src_io
  ../../src_io/numerical_basis.cpp ../../src_io/bessel_basis.cpp ../../src_io/numerical_descriptor.cpp ../../src_io/output.cpp
  ../../src_io/dos.cpp ../../src_io/energy_dos_pw.cpp ../../src_io/wf_io.cpp ../../src_io/rwstream.cpp
  ../../src_io/eximport.cpp ../../src_io/write_input.cpp ../../src_io/read_rho.cpp ../../src_io/write_rho.cpp
  ../../src_io/write_pot.cpp  ../../src_io/write_rho_cube.cpp  ../../src_io/write_rho_dipole.cpp ../../src_io/optical.cpp
  ../../src_io/cal_test.cpp  ../../src_io/winput.cpp ../../src_io/chi0_hilbert.cpp
  ../../src_io/chi0_standard.cpp  ../../src_io/epsilon0_pwscf.cpp  ../../src_io/epsilon0_vasp.cpp
  ../../src_io/to_wannier90.cpp  ../../src_io/unk_overlap_pw.cpp ../../src_io/berryphase.cpp
  ../../src_io/print_info.cpp ../../src_io/restart.cpp
  # module_cell
  ../../module_cell/atom_spec.cpp ../../module_cell/unitcell.cpp   ../../module_cell/read_pp.cpp ../../module_cell/read_atoms.cpp  ../../module_cell/read_cell_pseudopots.cpp
  ../../module_cell/read_pp_upf100.cpp  ../../module_cell/read_pp_upf201.cpp  ../../module_cell/read_pp_vwr.cpp  ../../module_cell/read_pp_blps.cpp
  ../../module_cell/atom_pseudo.cpp  ../../module_cell/pseudo_nc.cpp ../../module_cell/unitcell_pseudo.cpp ../../module_cell/setup_nonlocal.cpp
  # module_md #In fact, md module should not contain any LCAO information. And then we can use lib md directly.
  ../../module_md/run_md_classic.cpp ../../module_md/LJ_potential.cpp ../../module_md/DP_potential.cpp ../../module_md/MD_func.cpp
  ../../module_md/cmd_neighbor.cpp ../../module_md/verlet.cpp ../../module_md/NVE.cpp ../../module_md/MSST.cpp 
  ../../module_md/NVT_ADS.cpp ../../module_md/NVT_NHC.cpp ../../module_md/FIRE.cpp ../../module_md/Langevin.cpp
  # module_re #In fact, relax module should not contain any LCAO information. And then we can use lib relax directly.
  ../../module_relaxation/bfgs_basic.cpp  ../../module_relaxation/ions.cpp  ../../module_relaxation/ions_move_basic.cpp  ../../module_relaxation/ions_move_bfgs.cpp
  ../../module_relaxation/ions_move_cg.cpp  ../../module_relaxation/ions_move_methods.cpp  ../../module_relaxation/ions_move_sd.cpp
  ../../module_relaxation/lattice_change_basic.cpp  ../../module_relaxation/lattice_change_cg.cpp  ../../module_relaxation/lattice_change_methods.cpp
  ../../module_relaxation/variable_cell.cpp  ../../module_relaxation/relaxation.cpp
  # src_lcao
  ../../src_lcao/exx_lip.cpp #maybe it should be moved into module_xc
)

if(${ENABLE_LIBXC})
    target_link_libraries(abacus_pw.exe Libxc::xc)
endif()
