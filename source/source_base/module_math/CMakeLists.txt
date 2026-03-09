# Module: Math
# Mathematical functions and utilities

set(MODULE_NAME module_math)

set(HEADERS
    math_bspline.h
    math_chebyshev.h
    math_erf_complex.h
    math_integral.h
    math_lebedev_laikov.h
    math_polyint.h
    math_sphbes.h
    math_ylmreal.h
    mymath.h
)

set(SOURCES
    math_bspline.cpp
    math_chebyshev.cpp
    math_erf_complex.cpp
    math_integral.cpp
    math_lebedev_laikov.cpp
    math_polyint.cpp
    math_sphbes.cpp
    math_ylmreal.cpp
    mymath.cpp
)

add_library(${MODULE_NAME} OBJECT ${SOURCES} ${HEADERS})

target_include_directories(${MODULE_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/source/source_base
)

target_link_libraries(${MODULE_NAME}
    PUBLIC
        base
        container
        ${math_libs}
        ${FFTW_LIBRARIES}
)

if(BUILD_TESTING)
    add_subdirectory(test)
endif()