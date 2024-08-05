list(APPEND objects
    cal_ddpsir_ylm.cpp
    cal_dpsir_ylm.cpp
    cal_dpsirr_ylm.cpp
    cal_psir_ylm.cpp
    gint_tools.cpp
    init_orb.cpp
    gint_inout.hpp)

add_library(
    tools
    SHARED
    ${objects}
)