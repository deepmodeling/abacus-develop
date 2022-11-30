list(APPEND objects
    berryphase.cpp
    bessel_basis.cpp
    cal_test.cpp
    dos.cpp
    energy_dos_pw.cpp
    numerical_basis.cpp
    numerical_descriptor.cpp
    output.cpp
    print_info.cpp
    read_rho.cpp
    restart.cpp
    rwstream.cpp
    to_wannier90.cpp
    unk_overlap_pw.cpp
    wf_io.cpp
    winput.cpp
    write_input.cpp    
    write_rho.cpp
    write_rho_cube.cpp
    write_wfc_realspace.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      cal_r_overlap_R.cpp
      energy_dos.cpp
      istate_charge.cpp
      istate_envelope.cpp
      mulliken_charge.cpp
      read_dm.cpp
      unk_overlap_lcao.cpp
      wf_local.cpp
      write_HS.cpp
      write_HS_R.cpp
      write_dm.cpp
  )
endif()

add_library(
    io
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(io)
endif()
