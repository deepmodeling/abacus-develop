remove_definitions(-D__MPI)

install(FILES sphbes_ref.dat func.dat qlist.dat fq_ref.dat
  qlist1.dat fq.dat d2fq_ref.dat qlist2.dat fq_fit_ref.dat
  H.LDA_PW-JTH.xml Fe.GGA_PBE-JTH.xml O.GGA_PBE-JTH.xml
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET Test_Read_Paw
  LIBS ${math_libs} base device
  SOURCES test_paw.cpp ../paw_element.cpp ../paw_sphbes.cpp ../paw_interface.cpp
)

AddTest(
  TARGET Test_SphBes_Func
  LIBS ${math_libs} base device
  SOURCES test_paw.cpp ../paw_element.cpp ../paw_sphbes.cpp ../paw_interface.cpp
)

AddTest(
  TARGET Test_SphBes_Transform
  LIBS ${math_libs} base device
  SOURCES test_paw.cpp ../paw_element.cpp ../paw_sphbes.cpp ../paw_interface.cpp
)

AddTest(
  TARGET Test_Spline
  LIBS ${math_libs} base device
  SOURCES test_paw.cpp ../paw_element.cpp ../paw_sphbes.cpp ../paw_interface.cpp
)

AddTest(
  TARGET Test_Paw_Cell
  LIBS ${math_libs} base device
  SOURCES test_paw1.cpp ../paw_element.cpp ../paw_sphbes.cpp ../paw_interface.cpp ../paw_cell.cpp
)