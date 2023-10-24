add_subdirectory(operator_pw)

list(APPEND objects
    hamilt_pw.cpp
    operator_pw/ekinetic_pw.cpp
    operator_pw/veff_pw.cpp
    operator_pw/nonlocal_pw.cpp
    operator_pw/meta_pw.cpp
    operator_pw/velocity_pw.cpp
    operator_pw/operator_pw.cpp
    forces.cpp
    stress_func_cc.cpp
    stress_func_ewa.cpp
    stress_func_gga.cpp
    stress_func_mgga.cpp
    stress_func_har.cpp
    stress_func_kin.cpp
    stress_func_loc.cpp
    stress_func_nl.cpp
    stress_pw.cpp
    VL_in_pw.cpp
    VNL_in_pw.cpp
    VNL_grad_pw.cpp
    wavefunc.cpp
    wf_atomic.cpp
    structure_factor.cpp
    structure_factor_k.cpp
    soc.cpp
    global.cpp
    parallel_grid.cpp
)

add_library(
    hamilt_pwdft
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hamilt_pwdft)
endif()

if (BUILD_TESTING)
    add_subdirectory(test)
    add_subdirectory(kernels/test)
endif()
