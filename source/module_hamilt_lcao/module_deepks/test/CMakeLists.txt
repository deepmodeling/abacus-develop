remove_definitions(-DUSE_PAW)

add_executable(
    test_deepks
    main_deepks.cpp klist_1.cpp LCAO_deepks_test_prep.cpp LCAO_deepks_test.cpp nnr.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(test_deepks)
endif()

get_target_property(ABACUS_LINK_LIBRARIES ${ABACUS_BIN_NAME} LINK_LIBRARIES)
target_link_libraries(
    test_deepks
    base cell symmetry md surchem xc_
    neighbor orb io_basic io_advanced relax gint driver esolver hsolver psi elecstate
    hamilt_general hamilt_pwdft hamilt_lcao tddft hamilt_ofdft hamilt_stodft planewave
    pthread vdw dftu sc hcontainer
    deepks device numerical_atomic_orbitals container psi_initializer
    ${ABACUS_LINK_LIBRARIES}
)
if(USE_ELPA)
  target_link_libraries(
      test_deepks
      genelpa
  )
endif()
if(USE_CUSOLVER_LCAO)
  target_link_libraries(diag_cusolver)
endif()

if (ENABLE_LIBRI)
  target_link_libraries(test_deepks
          ri)
endif()

install(
    TARGETS test_deepks
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../../../tests
)
