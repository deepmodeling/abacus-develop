list(APPEND objects
    H_Ewald_pw.cpp
    VL_in_pw.cpp
    VNL_in_pw.cpp
    VNL_grad_pw.cpp
    charge.cpp
    charge_broyden.cpp
    charge_extra.cpp
    charge_mixing.cpp
    charge_pulay.cpp
    energy.cpp
    forces.cpp
    global.cpp
    klist.cpp
    magnetism.cpp
    occupy.cpp
    structure_factor.cpp
    pw_complement.cpp
    soc.cpp
    sto_iter.cpp
    sto_hchi.cpp
    sto_wf.cpp
    sto_func.cpp
    sto_forces.cpp
    sto_stress_pw.cpp
    stress_func_cc.cpp
    stress_func_ewa.cpp
    stress_func_gga.cpp
    stress_func_mgga.cpp
    stress_func_har.cpp
    stress_func_kin.cpp
    stress_func_loc.cpp
    stress_func_nl.cpp
    stress_func_print.cpp
    stress_pw.cpp
    of_stress_pw.cpp
    symmetry_rho.cpp
    threshold_elec.cpp
    wavefunc.cpp
    wf_atomic.cpp
    wf_igk.cpp
    xc_3.cpp
)

add_library(
    pw
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(pw)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

