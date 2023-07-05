AddTest(
  TARGET numerical_radial
  SOURCES
    numerical_radial_test.cpp
    ../numerical_radial.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET atomic_radials
  SOURCES
    atomic_radials_test.cpp
    ../atomic_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET beta_radials
  SOURCES
    beta_radials_test.cpp
    ../beta_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET radial_collection
  SOURCES
    radial_collection_test.cpp
    ../radial_collection.cpp
    ../atomic_radials.cpp
    ../beta_radials.cpp
    ../radial_set.cpp
    ../numerical_radial.cpp
  LIBS ${math_libs} device base
)

AddTest(
  TARGET real_gaunt_table
  SOURCES
    real_gaunt_table_test.cpp
    ../real_gaunt_table.cpp
    ../../module_ao/ORB_gaunt_table.cpp
  LIBS ${math_libs} device base container kernels
)

