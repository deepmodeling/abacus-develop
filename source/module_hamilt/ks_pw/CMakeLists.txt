list(APPEND operator_ks_pw_srcs
    operator_pw.cpp
    ekinetic_pw.cpp
    veff_pw.cpp
    nonlocal_pw.cpp
    meta_pw.cpp
    velocity_pw.cpp
    ../operator.cpp
    ../src/ekinetic.cpp
    ../src/nonlocal.cpp
    ../src/veff.cpp
)

if (USE_CUDA) 
    list(APPEND operator_ks_pw_srcs ../src/cuda/ekinetic.cu ../src/cuda/nonlocal.cu ../src/cuda/veff.cu)
elseif(USE_ROCM)
    list(APPEND operator_ks_pw_srcs ../src/rocm/ekinetic.cu ../src/rocm/nonlocal.cu ../src/rocm/veff.cu)
endif()

add_library(
    operator_ks_pw
    OBJECT
    ${operator_ks_pw_srcs}
)

if(ENABLE_COVERAGE)
  add_coverage(operator_ks_pw)
endif()
