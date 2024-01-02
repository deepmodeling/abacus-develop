AddTest(
  TARGET numerical_radial
  SOURCES
    numerical_radial_test.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET atomic_radials
  SOURCES
    atomic_radials_test.cpp
    ../atomic_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../../module_io/output_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET beta_radials
  SOURCES
    beta_radials_test.cpp
    ../beta_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../../module_io/output_radial.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET hydrogen_radials
  SOURCES
    hydrogen_radials_test.cpp
    ../hydrogen_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../../module_io/output_radial.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET pswfc_radials
  SOURCES
  pswfc_radials_test.cpp
    ../pswfc_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
    ../../../module_io/output_radial.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET sphbes_radials
  SOURCES
  sphbes_radials_test.cpp
    ../sphbes_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../../module_io/output_radial.cpp
    ../../module_ao/ORB_atomic_lm.cpp
 LIBS ${math_libs} device base
)

AddTest(
  TARGET radial_collection
  SOURCES
    radial_collection_test.cpp
    ../radial_collection.cpp
    ../atomic_radials.cpp
    ../beta_radials.cpp
    ../hydrogen_radials.cpp
    ../pswfc_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
    ../../../module_io/output_radial.cpp
    ../sphbes_radials.cpp
    ../../module_ao/ORB_atomic_lm.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET two_center_table
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
    ../../../module_io/output_radial.cpp
  LIBS ${math_libs} device base container orb
)

AddTest(
  TARGET real_gaunt_table
  SOURCES
    real_gaunt_table_test.cpp
    ../real_gaunt_table.cpp
    ../../module_ao/ORB_gaunt_table.cpp
  LIBS ${math_libs} device base container
)

AddTest(
  TARGET two_center_integrator
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
    ../../../module_io/output_radial.cpp
  LIBS ${math_libs} device base container orb
)

AddTest(
  TARGET two_center_bundle
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
    ../../../module_io/output_radial.cpp
  LIBS ${math_libs} device base container orb
)

