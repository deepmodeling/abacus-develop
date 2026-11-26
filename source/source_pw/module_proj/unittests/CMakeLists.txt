addtest(
  TARGET MODULE_PW_radial_proj_test
  LIBS parameter base device
  SOURCES test_radial_proj.cpp ../radial_proj.cpp)

addtest(
  TARGET MODULE_PW_projgen_test
  LIBS parameter base device
  SOURCES test_projgen.cpp ../projgen.cpp)

addtest(
  TARGET MODULE_PW_onsite_proj_init_test
  LIBS parameter
       base
       device
       psi
       planewave
       cell
       domain_decomposition
       neighbor
       neighbor_search
       symmetry
       nao_smoothgen
  SOURCES test_onsite_proj_init.cpp ../onsite_proj.cpp ../onsite_proj_init.cpp
          ../onsite_proj_print.cpp ../radial_proj.cpp)
set_tests_properties(MODULE_PW_onsite_proj_init_test
                     PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
