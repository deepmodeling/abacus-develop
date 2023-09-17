list(APPEND objects
    gint.cpp
    gint_gamma.cpp
    gint_gamma_env.cpp
    gint_gamma_vl.cpp
    gint_fvl.cpp
    gint_rho.cpp
    gint_tau.cpp
    gint_vl.cpp
    gint_k_env.cpp
    gint_k_sparse.cpp
    gint_k_sparse1.cpp
    gint_k_pvpr.cpp
    gint_k_pvdpr.cpp
    gint_tools.cpp
    grid_bigcell.cpp
    grid_meshball.cpp
    grid_meshcell.cpp
    grid_meshk.cpp
    grid_technique.cpp
    gint_vl.cu
)

add_library(
    gint
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(gint)
endif()
