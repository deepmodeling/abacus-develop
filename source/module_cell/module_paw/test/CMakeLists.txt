remove_definitions(-D__MPI)

install(FILES sphbes_ref.dat func.dat qlist.dat fq_ref.dat
  qlist1.dat fq.dat d2fq_ref.dat qlist2.dat fq_fit_ref.dat
  kpg.dat ylm_ref.dat leg_pol.dat kpg1.dat igxyz.dat
  eigts.dat ca.dat rhoij.dat vkb_ref.dat
  Si_test.xml gnorm.dat ptilde_ref.dat psi.dat rhoij1.dat vnlpsi_ref.dat
  H.LDA_PW-JTH.xml Fe.GGA_PBE-JTH.xml O.GGA_PBE-JTH.xml Si.GGA_PBE-JTH.xml
  STRU snlpsi_ref.dat
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_PAW)
  install(FILES C.xml H.xml 
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

AddTest(
  TARGET Test_Paw
  LIBS ${math_libs} base device
  SOURCES test_paw.cpp ../paw_element.cpp ../paw_sphbes.cpp ../paw_atom.cpp
)

AddTest(
  TARGET Test_Paw1
  LIBS ${math_libs} base device
  SOURCES test_paw1.cpp ../paw_element.cpp ../paw_sphbes.cpp ../paw_cell.cpp ../paw_atom.cpp
)

AddTest(
  TARGET Test_Paw2
  LIBS ${math_libs} base device
  SOURCES test_paw2.cpp ../paw_atom.cpp
)

if(ENABLE_PAW)
  AddTest(
    TARGET Test_Paw3
    LIBS ${math_libs} base device libpaw_interface -lifcore
    SOURCES test_paw3.cpp ../paw_atom.cpp ../paw_cell_libpaw.cpp
  )
  add_dependencies(Test_Paw3 libpaw_interface)
else()
  AddTest(
    TARGET Test_Paw3
    LIBS ${math_libs} base device
    SOURCES test_paw3.cpp ../paw_atom.cpp ../paw_cell_libpaw.cpp
  )
endif()