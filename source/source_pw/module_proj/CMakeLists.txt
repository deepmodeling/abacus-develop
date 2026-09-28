# module_proj: projection of plane-wave wavefunctions onto local (atomic) orbitals.
#
# This module hosts the onsite-projector code extracted from module_pwdft:
# radial projectors (SBT), onsite projector init/tabulate/overlap, and the
# Onsite_Proj_tools computation backend.
#
# Known temporary dependency (recorded during extraction):
# - Onsite_Proj_tools includes source_pw/module_pwdft/vnl_pw.h and
#   source_pw/module_pwdft/kernels/{force_op,stress_op}.h, i.e. module_proj
#   currently depends on module_pwdft headers. This is tolerated short-term;
#   the planned cleanup is to split the pseudopotential-specific constructor
#   of Onsite_Proj_tools so this reverse dependency can be removed.

list(APPEND objects
    radial_proj.cpp
    onsite_proj_print.cpp
    onsite_proj.cpp
    onsite_proj_init.cpp
    onsite_proj_overlap.cpp
    onsite_proj_tools.cpp
    onsite_proj_tools_becp.cpp
    onsite_proj_tools_dbecp.cpp
    onsite_proj_tools_force.cpp
    onsite_proj_tools_stress.cpp
)

add_library(
    module_proj
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(module_proj)
endif()

if(BUILD_TESTING)
  add_subdirectory(unittests)
endif()
