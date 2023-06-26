remove_definitions(-D__MPI)

install(FILES H.LDA_PW-JTH.xml sphbes_ref.dat func.dat qlist.dat fq_ref.dat
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET Test_Read_Paw
  LIBS ${math_libs} base device
  SOURCES test_paw.cpp ../paw_element.cpp ../paw_sphbes.cpp
)

AddTest(
  TARGET Test_SphBes_Func
  LIBS ${math_libs} base device
  SOURCES test_paw.cpp ../paw_element.cpp ../paw_sphbes.cpp
)

AddTest(
  TARGET Test_SphBes_Transform
  LIBS ${math_libs} base device
  SOURCES test_paw.cpp ../paw_element.cpp ../paw_sphbes.cpp
)