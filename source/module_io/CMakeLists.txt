list(APPEND objects
    input.cpp
    input_conv.cpp
    bessel_basis.cpp
    cal_test.cpp
    dos.cpp
    write_dos_pw.cpp
    nscf_band.cpp
    write_istate_info.cpp
    numerical_basis.cpp
    numerical_descriptor.cpp
    output.cpp
    print_info.cpp
    read_rho.cpp
    restart.cpp
    rwstream.cpp
    unk_overlap_pw.cpp
    write_wfc_pw.cpp
    write_input.cpp
    write_rho.cpp
    write_dipole.cpp
    write_wfc_r.cpp
)

list(APPEND objects_advanced
    berryphase.cpp
    to_wannier90.cpp
    winput.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      cal_r_overlap_R.cpp
      write_dos_lcao.cpp
      write_orb_info.cpp
      write_proj_band_lcao.cpp
      nscf_fermi_surf.cpp
      istate_charge.cpp
      istate_envelope.cpp
      read_dm.cpp
      unk_overlap_lcao.cpp
      read_wfc_nao.cpp
      write_wfc_nao.cpp
      write_HS.cpp
      write_HS_sparse.cpp
      write_HS_R.cpp
      write_dm.cpp
      write_dm_sparse.cpp
  )
  list(APPEND objects_advanced
      mulliken_charge.cpp
  )
endif()

add_library(
    io
    OBJECT
    ${objects} ${objects_advanced}
)

add_library(
    io_cov
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(io_cov)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
