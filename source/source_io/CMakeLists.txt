# Note: Help data is embedded in Input_Item registrations under module_parameter.
# The old input_help_data.h generation has been removed.

list(
  APPEND
  objects
  module_parameter/input_conv.cpp
  module_ctrl/ctrl_output_fp.cpp
  module_ctrl/ctrl_output_pw.cpp
  module_ctrl/ctrl_output_td.cpp
  module_ctrl/ctrl_output_td_pw.cpp
  module_bessel/bessel_basis.cpp
  module_output/cal_test.cpp
  module_dos/cal_dos.cpp
  module_dos/cal_ldos.cpp
  module_ml/cal_mlkedf_desc.cpp
  module_dos/write_dos_pw.cpp
  module_energy/write_bands.cpp
  module_energy/write_eig_occ.cpp
  module_energy/td_eij_io.cpp
  module_energy/td_eij_pw.cpp
  module_bessel/numerical_basis.cpp
  module_bessel/numerical_basis_overlap.cpp
  module_bessel/numerical_basis_output.cpp
  module_bessel/numerical_basis_jyjy.cpp
  module_bessel/numerical_descriptor.cpp
  module_output/print_info.cpp
  module_output/read_cube.cpp
  module_wf/read_wfc_pw.cpp
  module_wf/read_wf2rho_pw.cpp
  module_restart/restart.cpp
  module_wf/write_wfc_pw.cpp
  module_output/write_pao.cpp
  module_output/write_cube.cpp
  module_elf/write_elf.cpp
  module_dipole/write_dipole.cpp
  module_ml/write_mlkedf_desc.cpp
  module_current/td_current_io.cpp
  module_current/td_current_io_comm.cpp
  module_current/td_current_pw.cpp
  module_current/td_spin_pw_io.cpp
  module_efield/td_efield_io.cpp
  module_efield/td_vecpot_io.cpp
  module_chgpot/write_libxc_r.cpp
  module_output/output_log.cpp
  module_output/td_conserv_io.cpp
  module_json/para_json.cpp
  parse_args.cpp
  input_help.cpp
  module_chgpot/get_pchg_pw.cpp
  module_wf/get_wf_pw.cpp)

list(
  APPEND
  objects_advanced
  module_unk/unk_overlap_pw.cpp
  module_unk/berryphase.cpp
  module_wannier/to_w90.cpp
  module_wannier/to_w90_pw.cpp
  module_wannier/to_w90_pw_io.cpp
  module_wannier/to_w90_pw_setup.cpp
  module_wannier/to_w90_pw_overlap.cpp
  module_wannier/to_w90_lcao_pw.cpp
  module_wannier/to_w90_lcao.cpp
  module_wannier/to_w90_lcao_io.cpp
  module_wannier/to_w90_lcao_setup.cpp
  module_wannier/to_w90_lcao_overlap.cpp
  module_wannier/fr_overlap.cpp)

if(ENABLE_LCAO)
  list(
    APPEND
    objects
    module_dos/write_dos_lcao.cpp
    module_dos/cal_pdos_gamma.cpp
    module_dos/cal_pdos_multik.cpp
    module_dos/write_pdos_text.cpp
    module_output/write_orb_info.cpp
    module_energy/write_proj_band_lcao.cpp
    module_chgpot/get_pchg_lcao.cpp
    module_wf/get_wf_lcao.cpp
    module_wf/read_wfc_nao.cpp
    module_wf/write_wfc_nao.cpp
    module_dm/write_dmk.cpp
    module_dm/write_dmr.cpp
    module_qo/to_qo_kernel.cpp
    module_qo/to_qo_mpi.cpp
    module_qo/to_qo_structures.cpp
    module_mulliken/output_sk.cpp
    module_mulliken/output_dmk.cpp
    module_mulliken/output_mulliken.cpp
    module_ml/io_npz.cpp
    module_hs/angmom_op_mat.cpp
    module_dhs/write_dh.cpp
    module_dhs/write_dh_terms.cpp)
  list(
    APPEND
    objects_advanced
    module_unk/unk_overlap_lcao.cpp
    module_hs/hs_dense_io.cpp
    module_hs/hs_sparse_io.cpp
    module_hs/dhs_sparse_writer.cpp
    module_hs/hsk_writer.cpp
    module_hs/hsr_writer.cpp
    module_hs/hsr_legacy.cpp
    module_hs/hterm_writer.cpp
    module_hs/lat_r_csr.cpp
    module_hs/pos_op_csr.cpp
    module_hs/pos_op_basis.cpp
    module_hs/pos_op_calc.cpp
    module_hs/pos_op_writer.cpp
    module_hs/pos_op_mat.cpp
    module_hs/vxc_op_tools.cpp
    module_hs/vxc_op_mat.cpp
    module_hs/vxc_op_r.cpp
    module_hs/vxc_op_lip.cpp
    module_hs/output_mat_sparse.cpp
    module_energy/write_eband_terms.cpp
    module_ctrl/ctrl_scf_lcao.cpp
    module_ctrl/ctrl_runner_lcao.cpp
    module_ctrl/ctrl_iter_lcao.cpp)
  if(ENABLE_LIBRI)
    list(APPEND objects module_restart/restart_exx_csr.cpp)
  endif()
endif()

add_library(
  io_input OBJECT
  module_parameter/read_inp_sys.cpp
  module_parameter/read_inp_estruc.cpp
  module_parameter/read_inp_relax.cpp
  module_parameter/read_input_item_md.cpp
  module_parameter/read_inp_ofdft.cpp
  module_parameter/read_input_item_sdft.cpp
  module_parameter/read_inp_tddft.cpp
  module_parameter/read_inp_bse.cpp
  module_parameter/read_inp_dfpt.cpp
  module_parameter/read_inp_deepks.cpp
  module_parameter/read_inp_model.cpp
  module_parameter/read_inp_sccs.cpp
  module_parameter/sccs_inp_check.cpp
  module_parameter/read_inp_postproc.cpp
  module_parameter/read_inp_exx_dftu.cpp
  module_parameter/read_inp_other.cpp
  module_parameter/read_inp_out.cpp
  module_parameter/availability.cpp
  module_parameter/availability_validator.cpp
  module_parameter/input_domain.cpp
  module_parameter/read_input.cpp
  module_parameter/read_set_globalv.cpp)

# Prints docs/parameters.yaml for the documentation build. It links only the
# INPUT registry and the few source_base files it calls, so .readthedocs.yaml
# builds this target instead of the full executable.
add_executable(
  abacus_input_yaml EXCLUDE_FROM_ALL
  input_yaml_main.cpp
  input_help.cpp
  ../source_base/global_file.cpp
  ../source_base/global_function.cpp
  ../source_base/global_variable.cpp
  ../source_base/parallel_common.cpp
  ../source_base/parallel_comm.cpp
  ../source_base/parallel_global.cpp
  ../source_base/parallel_reduce.cpp
  ../source_base/timer.cpp
  ../source_base/tool_quit.cpp
  ../source_base/tool_title.cpp
  ../source_base/module_device/device.cpp)
target_link_libraries(abacus_input_yaml PRIVATE io_input parameter
                                                abacus::linalg_libs)

add_library(io_basic OBJECT ${objects})

# Include directory for build_info.h (generated at build time)
target_include_directories(io_basic PUBLIC ${CMAKE_BINARY_DIR}/source/source_io)

add_library(io_advanced OBJECT ${objects_advanced})

if(ENABLE_COVERAGE)
  add_coverage(io_basic)
endif()

if(BUILD_TESTING)
  add_subdirectory(module_parameter/test)
  add_subdirectory(module_hs/unittests)
  add_subdirectory(test_serial)
  add_subdirectory(unittests)
  add_subdirectory(module_dos/unittests)
endif()

if(ENABLE_JSON)
  add_subdirectory(module_json)
endif()
