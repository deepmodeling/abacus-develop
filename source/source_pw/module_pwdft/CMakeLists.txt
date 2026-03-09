

list(APPEND objects
    hamilt_pw.cpp
    op_pw_ekin.cpp
    op_pw_veff.cpp
    op_pw_nl.cpp
    op_pw_meta.cpp
    op_pw_vel.cpp
    op_pw.cpp
    op_pw_proj.cpp
    op_pw_exx.cpp
    op_pw_exx_ace.cpp
    op_pw_exx_pot.cpp
    setup_pot.cpp
    setup_pwrho.cpp
    setup_pwwfc.cpp
    update_cell_pw.cpp
    dftu_pw.cpp
    deltaspin_pw.cpp
    forces_nl.cpp
    forces_cc.cpp
    forces_scc.cpp
    forces.cpp
    forces_us.cpp
    forces_onsite.cpp
    stress_cc.cpp
    stress_ewa.cpp
    stress_gga.cpp
    stress_mgga.cpp
    stress_har.cpp
    stress_kin.cpp
    stress_loc.cpp
    stress_nl.cpp
    stress_us.cpp
    stress_onsite.cpp
    stress_exx.cpp
    stress_pw.cpp
    vl_pw.cpp
    vnl_pw.cpp
    vnl_pw_grad.cpp
    structure_factor.cpp
    structure_factor_k.cpp
    soc.cpp
    parallel_grid.cpp
    elecond.cpp
    fs_nonlocal_tools.cpp
    fs_kin_tools.cpp
    radial_proj.cpp
    onsite_projector.cpp
    onsite_proj_tools.cpp
    vsep_pw.cpp
)

add_library(
    module_pwdft
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(module_pwdft)
endif()

if (BUILD_TESTING)
    if(ENABLE_MPI)
        add_subdirectory(test)
        add_subdirectory(kernels/test)
    endif()
endif()
