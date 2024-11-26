list(APPEND objects
    input_conv.cpp
    bessel_basis.cpp
    cal_test.cpp
    cal_dos.cpp
    cif_io.cpp
    write_dos_pw.cpp
    nscf_band.cpp
    write_istate_info.cpp
    numerical_basis.cpp
    numerical_basis_jyjy.cpp
    numerical_descriptor.cpp
    output.cpp
    print_info.cpp
    read_cube.cpp
    rhog_io.cpp
    read_exit_file.cpp
    read_wfc_pw.cpp
    read_wfc_to_rho.cpp
    restart.cpp
    binstream.cpp
    write_wfc_pw.cpp
    write_cube.cpp
    write_elecstat_pot.cpp
    write_elf.cpp
    write_dipole.cpp
    td_current_io.cpp
    write_wfc_r.cpp
    output_log.cpp
    para_json.cpp
    parse_args.cpp
    orb_io.cpp
)

list(APPEND objects_advanced
    unk_overlap_pw.cpp
    berryphase.cpp
    to_wannier90.cpp
    to_wannier90_pw.cpp
    to_wannier90_lcao_in_pw.cpp
    to_wannier90_lcao.cpp
    fR_overlap.cpp
    winput.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      write_dos_lcao.cpp
      write_orb_info.cpp
      write_proj_band_lcao.cpp
      nscf_fermi_surf.cpp
      get_pchg_lcao.cpp
      get_wf_lcao.cpp
      read_wfc_nao.cpp
      read_wfc_lcao.cpp
      write_wfc_nao.cpp
      io_dmk.cpp
      write_dmr.cpp
      dos_nao.cpp
      sparse_matrix.cpp
      file_reader.cpp
      csr_reader.cpp
      to_qo_kernel.cpp
      to_qo_mpi.cpp
      to_qo_structures.cpp
      output_sk.cpp
      output_dmk.cpp
      output_mulliken.cpp
      io_npz.cpp
  )
  list(APPEND objects_advanced
      unk_overlap_lcao.cpp
      write_HS_R.cpp
      write_HS_sparse.cpp
      single_R_io.cpp
      cal_r_overlap_R.cpp
      output_mat_sparse.cpp
  )
endif()

add_library(
    io_input
    OBJECT
    read_input_item_system.cpp
    read_input_item_elec_stru.cpp
    read_input_item_relax.cpp
    read_input_item_md.cpp
    read_input_item_ofdft.cpp
    read_input_item_sdft.cpp
    read_input_item_tddft.cpp
    read_input_item_deepks.cpp
    read_input_item_model.cpp
    read_input_item_postprocess.cpp  
    read_input_item_exx_dftu.cpp
    read_input_item_other.cpp   
    read_input_item_output.cpp  
    read_input.cpp
    read_set_globalv.cpp
)

add_library(
    io_basic
    OBJECT
    ${objects}
)

add_library(
    io_advanced
    OBJECT
    ${objects_advanced}
)

if(ENABLE_COVERAGE)
  add_coverage(io_basic)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
    add_subdirectory(test_serial)
  endif()
endif()

if(ENABLE_RAPIDJSON)
  if(ENABLE_MPI)
    add_subdirectory(json_output)
  endif()
endif()
