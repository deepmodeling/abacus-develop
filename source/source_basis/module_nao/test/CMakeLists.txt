AddTest(
  TARGET MODULE_NAO_numerical_radial
  SOURCES
    numerical_radial_test.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
  LIBS parameter ${math_libs} device base 
)

AddTest(
  TARGET MODULE_NAO_atomic_radials
  SOURCES
    atomic_radials_test.cpp
    ../atomic_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../module_ao/ORB_atomic.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base 
)

AddTest(
  TARGET MODULE_NAO_beta_radials
  SOURCES
    beta_radials_test.cpp
    ../beta_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../module_ao/ORB_atomic.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base 
)

AddTest(
  TARGET MODULE_NAO_hydrogen_radials
  SOURCES
    hydrogen_radials_test.cpp
    ../hydrogen_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../module_ao/ORB_atomic.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base 
)

AddTest(
  TARGET MODULE_NAO_pswfc_radials
  SOURCES
  pswfc_radials_test.cpp
    ../pswfc_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../module_ao/ORB_atomic.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base 
)

AddTest(
  TARGET MODULE_NAO_sphbes_radials
  SOURCES
  sphbes_radials_test.cpp
    ../sphbes_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../module_ao/ORB_atomic.cpp
    ../../../source_io/module_output/orb_io.cpp
 LIBS parameter ${math_libs} device base 
)

AddTest(
  TARGET MODULE_NAO_radial_collection
  SOURCES
    radial_collection_test.cpp
    ../radial_collection.cpp
    ../atomic_radials.cpp
    ../beta_radials.cpp
    ../hydrogen_radials.cpp
    ../pswfc_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../sphbes_radials.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../module_ao/ORB_atomic.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base 
)

AddTest(
  TARGET MODULE_NAO_two_center_table
  SOURCES
    two_center_table_test.cpp
    ../two_center_table.cpp
    ../radial_collection.cpp
    ../atomic_radials.cpp
    ../beta_radials.cpp
    ../hydrogen_radials.cpp
    ../pswfc_radials.cpp
    ../sphbes_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../two_center_bundle.cpp
    ../two_center_integrator.cpp
    ../real_gaunt_table.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base container orb 
)

AddTest(
  TARGET MODULE_NAO_real_gaunt_table
  SOURCES
    real_gaunt_table_test.cpp
    ../real_gaunt_table.cpp
    ../../module_ao/ORB_gaunt_table.cpp
  LIBS parameter ${math_libs} device base container
)

AddTest(
  TARGET MODULE_NAO_two_center_integrator
  SOURCES
    two_center_integrator_test.cpp
    ../two_center_integrator.cpp
    ../two_center_table.cpp
    ../real_gaunt_table.cpp
    ../radial_collection.cpp
    ../atomic_radials.cpp
    ../beta_radials.cpp
    ../hydrogen_radials.cpp
    ../pswfc_radials.cpp
    ../sphbes_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../two_center_bundle.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base container orb 
)

AddTest(
  TARGET MODULE_NAO_two_center_bundle
  SOURCES
    two_center_bundle_test.cpp
    ../two_center_bundle.cpp
    ../two_center_integrator.cpp
    ../two_center_table.cpp
    ../real_gaunt_table.cpp
    ../radial_collection.cpp
    ../atomic_radials.cpp
    ../beta_radials.cpp
    ../hydrogen_radials.cpp
    ../pswfc_radials.cpp
    ../sphbes_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base container orb 
)

AddTest(
  TARGET MODULE_NAO_generate_projection
  SOURCES
    projgen_test.cpp
    ../two_center_bundle.cpp
    ../two_center_integrator.cpp
    ../two_center_table.cpp
    ../real_gaunt_table.cpp
    ../radial_collection.cpp
    ../atomic_radials.cpp
    ../beta_radials.cpp
    ../hydrogen_radials.cpp
    ../pswfc_radials.cpp
    ../sphbes_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../../source_io/module_output/orb_io.cpp
  LIBS parameter ${math_libs} device base container orb 
)

