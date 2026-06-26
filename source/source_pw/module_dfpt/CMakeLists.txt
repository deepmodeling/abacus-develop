set(MODULE_NAME module_dfpt)

set(SOURCES
    dfpt_pw.cpp
    dfpt_pw_data.cpp
    dfpt_pert.cpp
    dfpt_stern.cpp
    dfpt_rho.cpp
    dfpt_phon.cpp
    dfpt_q0.cpp
    dfpt_metal.cpp
)

set(HEADERS
    dfpt_pw.h
    dfpt_pw_data.h
    dfpt_pert.h
    dfpt_stern.h
    dfpt_rho.h
    dfpt_phon.h
    dfpt_q0.h
    dfpt_metal.h
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