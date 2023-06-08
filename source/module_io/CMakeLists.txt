list(APPEND objects
    input.cpp
    input_conv.cpp
    bessel_basis.cpp
    cal_test.cpp
    cal_dos.cpp
    write_dos_pw.cpp
    nscf_band.cpp
    write_istate_info.cpp
    numerical_basis.cpp
    numerical_descriptor.cpp
    output.cpp
    print_info.cpp
    read_cube.cpp
    read_rho.cpp
    restart.cpp
    binstream.cpp
    write_wfc_pw.cpp
    write_input.cpp
    write_cube.cpp
    write_rho.cpp
    write_potential.cpp
    write_dipole.cpp
    write_wfc_r.cpp
    output_rho.cpp
    output_potential.cpp
)

list(APPEND objects_advanced
    unk_overlap_pw.cpp
    berryphase.cpp
    to_wannier90.cpp
    winput.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      write_dos_lcao.cpp
      write_orb_info.cpp
      write_proj_band_lcao.cpp
      nscf_fermi_surf.cpp
      istate_charge.cpp
      istate_envelope.cpp
      read_dm.cpp
      read_wfc_nao.cpp
      write_wfc_nao.cpp
      write_HS.cpp
      write_dm.cpp
      dos_nao.cpp
      output_dm.cpp
      output_dm1.cpp
  )
  list(APPEND objects_advanced
      unk_overlap_lcao.cpp
      mulliken_charge.cpp
      write_HS_R.cpp
      write_HS_sparse.cpp
      single_R_io.cpp
      write_dm_sparse.cpp
      cal_r_overlap_R.cpp
      output_mat_sparse.cpp
  )
endif()

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
  add_subdirectory(test)
  add_subdirectory(test_serial)
endif()
