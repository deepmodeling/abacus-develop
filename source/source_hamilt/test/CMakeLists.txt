AddTest(
  TARGET MODULE_HAMILT_ewald_dnrm2
  SOURCES dnrm2_test.cpp  ../module_ewald/dnrm2.cpp
)

AddTest(
  TARGET MODULE_HAMILT_ewald_rgen
  LIBS parameter ${math_libs} base device
  SOURCES rgen_test.cpp ../module_ewald/H_Ewald_pw.cpp ../module_ewald/dnrm2.cpp
)
