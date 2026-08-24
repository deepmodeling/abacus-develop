set(MODULE_NAME module_dfpt)

set(SOURCES
    dfpt_pw.cpp
    dfpt_pw_data.cpp
    dfpt_kq_basis.cpp
    dfpt_pert.cpp
    dfpt_stern.cpp
    dfpt_rho.cpp
    dfpt_phon.cpp
    dfpt_q0.cpp
    dfpt_metal.cpp
    dfpt_hamilt_shift.cpp
)

set(HEADERS
    dfpt_pw.h
    dfpt_pw_data.h
    dfpt_kq_basis.h
    dfpt_pert.h
    dfpt_stern.h
    dfpt_rho.h
    dfpt_phon.h
    dfpt_q0.h
    dfpt_metal.h
    dfpt_hamilt_shift.h
)

add_library(${MODULE_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${MODULE_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${MODULE_NAME}
    PUBLIC
        base
        cell
        psi
        elecstate
        module_pwdft
)

if (BUILD_TESTING)
    if(ENABLE_MPI)
        add_subdirectory(test)
        # serial (__MPI-off) tests for the FFT-driving kernels; must be built
        # inside ENABLE_MPI so the shared prebuilt test libs stay available
        add_subdirectory(test_serial)
    endif()
endif()